#!/usr/bin/env python3
"""
Run after qc_pipeline.py. Predicts 3D structures for a small shortlist of fully-valid heavy+light pairs
using IgFold, via a dedicated Python >=3.11 venv (IgFold is incompatible with
this skill's main venv's Python 3.9). The shortlist can be ranked by
abundance (default -- the most common pairs) or by fitness-LM perplexity
(--rank-by perplexity -- the pairs AbLang2/ESM-2 found least "natural",
useful for checking whether high-perplexity sequences actually show any
structural red flags, vs. perplexity just being a naturalness/novelty signal
unrelated to real foldability).

Structure prediction is far too slow to run at NGS scale, and genuinely
unnecessary for most reads -- this only ever scores a small, explicit
shortlist (default top 25 by abundance among reads passing both chains'
error threshold), never the full dataset or a large sample like the fitness
LM / KL divergence steps use. There's also no lightweight way to render a 3D
structure into a static PowerPoint image, so results land in the deck as a
table (id, abundance, mean predicted RMSD, PDB file path) via build_report.py
--structure-manifest, not as a picture.

One-time setup this depends on (see SKILL.md "Environment setup"):
    brew install python@3.11
    /opt/homebrew/bin/python3.11 -m venv venv-igfold311
    source venv-igfold311/bin/activate && pip install igfold

If that venv doesn't exist at --igfold-python, this script skips structure
prediction with a clear message rather than crashing -- it's an optional,
manually-bootstrapped prerequisite, not something every run can assume.

Example:
    python structure_prediction.py \\
        --scored-reads-csv antibody_qc_output/scored_reads.csv \\
        --heavy-col HAA --light-col LAA --read-count-col read_count \\
        --heavy-error-col Heavy_True_Errors --light-error-col Light_True_Errors \\
        --output-dir antibody_qc_output
"""

import argparse
import json
import os
import subprocess
import sys

import pandas as pd

DEFAULT_IGFOLD_PYTHON = os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "venv-igfold311", "bin", "python"))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--scored-reads-csv", required=True)
    parser.add_argument("--heavy-col", required=True)
    parser.add_argument("--light-col", required=True)
    parser.add_argument("--read-count-col", default="read_count")
    parser.add_argument("--heavy-error-col", default="Heavy_True_Errors")
    parser.add_argument("--light-error-col", default="Light_True_Errors")
    parser.add_argument("--error-threshold", type=int, default=2,
                         help="Same meaning as qc_pipeline.py -- only reads passing this on "
                              "both chains are eligible for the shortlist.")
    parser.add_argument("--shortlist-size", type=int, default=25,
                         help="How many fully-valid pairs to fold.")
    parser.add_argument("--rank-by", choices=["abundance", "perplexity"], default="abundance",
                         help="'abundance' (default) takes the most-common pairs by read count. "
                              "'perplexity' takes the highest-perplexity (least model-natural) "
                              "pairs instead -- requires --scored-reads-csv to already have a "
                              "perplexity column (i.e. point it at ml_bias_fitness.py's "
                              "scored_reads_with_fitness*.csv output, not the base scored_reads.csv).")
    parser.add_argument("--perplexity-col", default="Paired_Fitness_Perplexity",
                         help="Column to rank by when --rank-by perplexity. Defaults to the joint "
                              "pairing-compatibility score; use Heavy_Fitness_Perplexity or "
                              "Light_Fitness_Perplexity instead to rank by one chain alone.")
    parser.add_argument("--igfold-python", default=DEFAULT_IGFOLD_PYTHON,
                         help=f"Path to the dedicated IgFold venv's python interpreter. "
                              f"Default: {DEFAULT_IGFOLD_PYTHON}")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    manifest_path = os.path.join(args.output_dir, "structure_manifest.json")

    if not os.path.exists(args.igfold_python):
        print(f"IgFold venv not found at {args.igfold_python} -- skipping structure prediction. "
              f"See SKILL.md 'Environment setup' to create it.")
        with open(manifest_path, "w") as f:
            json.dump({"skipped": True, "reason": "igfold venv not found", "rows": []}, f, indent=2)
        print(f"Structure manifest (skipped): {manifest_path}")
        return

    df = pd.read_csv(args.scored_reads_csv)
    valid_both = df[(df[args.heavy_error_col] <= args.error_threshold) &
                    (df[args.light_error_col] <= args.error_threshold)]

    if args.rank_by == "perplexity":
        if args.perplexity_col not in valid_both.columns:
            parser.error(f"--rank-by perplexity needs column '{args.perplexity_col}' in "
                         f"--scored-reads-csv -- point it at ml_bias_fitness.py's "
                         f"scored_reads_with_fitness*.csv output.")
        valid_both = valid_both.dropna(subset=[args.perplexity_col])
        shortlist = valid_both.sort_values(args.perplexity_col, ascending=False).head(args.shortlist_size).copy()
        rank_note = f"by highest {args.perplexity_col}"
    else:
        shortlist = valid_both.sort_values(args.read_count_col, ascending=False).head(args.shortlist_size).copy()
        rank_note = "by abundance"

    shortlist["id"] = [f"seq{idx}" for idx in range(len(shortlist))]
    shortlist = shortlist.rename(columns={args.heavy_col: "heavy_seq", args.light_col: "light_seq"})

    print(f"Shortlisted {len(shortlist)} of {len(valid_both):,} fully-valid pairs "
          f"({rank_note}) for structure prediction.")

    extra_cols = [args.perplexity_col] if args.rank_by == "perplexity" else []
    structures_dir = os.path.join(args.output_dir, "structures")
    shortlist_csv = os.path.join(args.output_dir, "structure_shortlist.csv")
    results_csv = os.path.join(args.output_dir, "structure_results.csv")
    shortlist[["id", "heavy_seq", "light_seq", args.read_count_col] + extra_cols].to_csv(shortlist_csv, index=False)

    worker_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "igfold_worker.py")
    cmd = [args.igfold_python, worker_path,
           "--input-csv", shortlist_csv, "--id-col", "id",
           "--heavy-col", "heavy_seq", "--light-col", "light_seq",
           "--output-dir", structures_dir, "--results-csv", results_csv]
    print(f"Running IgFold worker: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    results = pd.read_csv(results_csv).merge(
        shortlist[["id", args.read_count_col] + extra_cols], on="id", how="left")

    if args.rank_by == "perplexity":
        headers = ["Sequence ID", "Read Count", args.perplexity_col, "Mean Predicted RMSD (Å)", "PDB Path"]
        rows = [
            [r["id"], f"{r[args.read_count_col]:,}", f"{r[args.perplexity_col]:.3f}",
             f"{r['mean_predicted_rmsd']:.3f}", r["pdb_path"]]
            for _, r in results.iterrows()
        ]
    else:
        headers = ["Sequence ID", "Read Count", "Mean Predicted RMSD (Å)", "PDB Path"]
        rows = [
            [r["id"], f"{r[args.read_count_col]:,}", f"{r['mean_predicted_rmsd']:.3f}", r["pdb_path"]]
            for _, r in results.iterrows()
        ]

    with open(manifest_path, "w") as f:
        json.dump({
            "skipped": False,
            "rank_by": args.rank_by,
            "headers": headers,
            "rows": rows,
        }, f, indent=2)
    print(f"Structure manifest: {manifest_path}")


if __name__ == "__main__":
    main()
