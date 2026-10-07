#!/usr/bin/env python3
"""
Run after qc_pipeline.py. Scores a random representative sample of reads for
humanness/immunogenicity risk using Sapiens, a human-antibody language
model, via a dedicated Python 3.7/3.8 environment (Sapiens is incompatible
with this skill's main venv's Python 3.9 -- see its own README).

Sapiens is a transformer forward pass per sequence like the fitness-LM step,
so the same auto-sampling logic applies (default 15,000 rows, same cost
reasoning as AbLang2 -- see ml_bias_fitness.py). Also scores each design's
own wild-type sequence as a reference baseline, same idea as
compute_wt_reference_perplexity / compute_ablang2_wt_reference elsewhere in
this skill.

One-time setup this depends on (see SKILL.md "Environment setup") -- Sapiens
needs Python <=3.8, which fails to build from source on modern macOS/arm64
(confirmed: the interpreter itself segfaults). Workaround used here: a
prebuilt (not locally-compiled) CPython distribution:
    curl -L -o /tmp/cpython38.tar.gz \\
        https://github.com/astral-sh/python-build-standalone/releases/download/20231002/cpython-3.8.18+20231002-aarch64-apple-darwin-install_only.tar.gz
    mkdir venv-sapiens38 && tar -xzf /tmp/cpython38.tar.gz -C /tmp && cp -R /tmp/python/* venv-sapiens38/
    venv-sapiens38/bin/python3.8 -m pip install sapiens pandas

If that environment doesn't exist at --sapiens-python, this script skips
humanness scoring with a clear message rather than crashing -- it's an
optional, manually-bootstrapped prerequisite, not something every run can
assume is present.

Example:
    python humanness_scoring.py \\
        --scored-reads-csv antibody_qc_output/scored_reads.csv \\
        --design-xlsx AffMat_design.xlsx \\
        --heavy-sheet "Heavy Chain Design" --heavy-col HAA \\
        --light-sheet "Light Chain Design" --light-col LAA \\
        --output-dir antibody_qc_output
"""

import argparse
import json
import os
import subprocess
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pwm_extraction import extract_full_pwm

DEFAULT_SAPIENS_PYTHON = os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "venv-sapiens38", "bin", "python3.8"))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--scored-reads-csv", required=True)
    parser.add_argument("--design-xlsx", required=True)
    parser.add_argument("--heavy-sheet", default=None)
    parser.add_argument("--heavy-col", default=None)
    parser.add_argument("--light-sheet", default=None)
    parser.add_argument("--light-col", default=None)
    parser.add_argument("--sample-size", type=int, default=15_000,
                         help="Random sample size for humanness scoring (same cost reasoning "
                              "as --fitness-sample-size in ml_bias_fitness.py). Pass 0 to force "
                              "scoring every row regardless of size.")
    parser.add_argument("--sample-seed", type=int, default=42)
    parser.add_argument("--sapiens-python", default=DEFAULT_SAPIENS_PYTHON,
                         help=f"Path to the dedicated Sapiens environment's python3.8. "
                              f"Default: {DEFAULT_SAPIENS_PYTHON}")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    if not args.heavy_sheet and not args.light_sheet:
        parser.error("Provide at least --heavy-sheet/--heavy-col.")

    os.makedirs(args.output_dir, exist_ok=True)
    manifest_path = os.path.join(args.output_dir, "humanness_manifest.json")

    if not os.path.exists(args.sapiens_python):
        print(f"Sapiens environment not found at {args.sapiens_python} -- skipping humanness "
              f"scoring. See SKILL.md 'Environment setup' to create it.")
        with open(manifest_path, "w") as f:
            json.dump({"skipped": True, "reason": "sapiens environment not found"}, f, indent=2)
        print(f"Humanness manifest (skipped): {manifest_path}")
        return

    df = pd.read_csv(args.scored_reads_csv)
    sampled = args.sample_size and len(df) > args.sample_size
    sample_df = df.sample(n=args.sample_size, random_state=args.sample_seed) if sampled else df

    chains = []
    if args.heavy_sheet:
        rules = extract_full_pwm(args.design_xlsx, args.heavy_sheet)
        chains.append(dict(label="Heavy", side="H", seq_col=args.heavy_col,
                            window_col="Heavy_Best_Window", rules=rules))
    if args.light_sheet:
        rules = extract_full_pwm(args.design_xlsx, args.light_sheet)
        chains.append(dict(label="Light", side="L", seq_col=args.light_col,
                            window_col="Light_Best_Window", rules=rules))

    bullets_data = {}
    worker_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sapiens_worker.py")

    for chain in chains:
        sample_note = f" (random sample of {len(sample_df):,} of {len(df):,} rows)" if sampled else ""
        print(f"Scoring {chain['label']} chain sequences with Sapiens{sample_note}...")

        # Score the aligned window, not the raw read column -- see the matching
        # comment in ml_bias_fitness.py. The WT reference sequences below are
        # already in the design's (short, V-gene-only) truncation, so the sample
        # needs to match that truncation too for the comparison to be meaningful.
        score_col = chain["window_col"] if chain["window_col"] in sample_df.columns else chain["seq_col"]
        wt_seqs = [rules["wt_sequence"] for rules in chain["rules"].values()]
        input_df = pd.DataFrame({
            "sequence": list(sample_df[score_col]) + wt_seqs,
            "chain": [chain["side"]] * (len(sample_df) + len(wt_seqs)),
        })
        input_csv = os.path.join(args.output_dir, f"_sapiens_input_{chain['label'].lower()}.csv")
        output_csv = os.path.join(args.output_dir, f"_sapiens_output_{chain['label'].lower()}.csv")
        input_df.to_csv(input_csv, index=False)

        subprocess.run([args.sapiens_python, worker_path,
                         "--input-csv", input_csv, "--seq-col", "sequence", "--chain-col", "chain",
                         "--output-csv", output_csv], check=True)

        results = pd.read_csv(output_csv)
        sample_results = results.iloc[:len(sample_df)]
        wt_results = results.iloc[len(sample_df):]

        mean_perplexity = sample_results["humanness_perplexity"].mean()
        mean_match_pct = sample_results["humanness_match_pct"].mean()
        wt_mean_perplexity = wt_results["humanness_perplexity"].mean()
        wt_mean_match_pct = wt_results["humanness_match_pct"].mean()

        print(f"  {chain['label']}: mean humanness perplexity {mean_perplexity:.3f} "
              f"(WT baseline {wt_mean_perplexity:.3f}), mean top-1 match {mean_match_pct:.1f}% "
              f"(WT baseline {wt_mean_match_pct:.1f}%)")

        bullets_data[chain["label"]] = {
            "mean_humanness_perplexity": round(float(mean_perplexity), 3),
            "wt_baseline_perplexity": round(float(wt_mean_perplexity), 3),
            "mean_top1_match_pct": round(float(mean_match_pct), 2),
            "wt_baseline_top1_match_pct": round(float(wt_mean_match_pct), 2),
            "rows_scored": len(sample_df),
        }

        os.remove(input_csv)
        os.remove(output_csv)

    summary_path = os.path.join(args.output_dir, "humanness_summary.json")
    with open(summary_path, "w") as f:
        json.dump(bullets_data, f, indent=2)
    print(f"Saved {summary_path}")

    chart_path = os.path.join(args.output_dir, "humanness_summary.png")
    plot_humanness_summary(bullets_data, chart_path)
    print(f"Saved {chart_path}")

    with open(manifest_path, "w") as f:
        json.dump({
            "skipped": False,
            "summary_json": summary_path,
            "summary_image": chart_path,
            "sampled": sampled,
            "rows_scored": len(sample_df),
            "total_rows": len(df),
        }, f, indent=2)
    print(f"Humanness manifest: {manifest_path}")


def plot_humanness_summary(bullets_data, output_path):
    """Grouped bar chart: sample mean vs. WT baseline, for both humanness
    metrics, one group per chain in `bullets_data`."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    chains = list(bullets_data.keys())
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

    x = np.arange(len(chains))
    width = 0.35

    ax1.bar(x - width / 2, [bullets_data[c]["mean_humanness_perplexity"] for c in chains],
            width, label="Sample Mean")
    ax1.bar(x + width / 2, [bullets_data[c]["wt_baseline_perplexity"] for c in chains],
            width, label="WT Baseline")
    ax1.set_xticks(x)
    ax1.set_xticklabels(chains)
    ax1.set_ylabel("Sapiens Humanness Perplexity (Lower = More Human-Like)", fontsize=9)
    ax1.set_title("Humanness Perplexity", fontsize=10, fontweight="bold")
    ax1.legend(fontsize=8)

    ax2.bar(x - width / 2, [bullets_data[c]["mean_top1_match_pct"] for c in chains],
            width, label="Sample Mean")
    ax2.bar(x + width / 2, [bullets_data[c]["wt_baseline_top1_match_pct"] for c in chains],
            width, label="WT Baseline")
    ax2.set_xticks(x)
    ax2.set_xticklabels(chains)
    ax2.set_ylabel("% Positions Matching Top-1 Human Residue", fontsize=9)
    ax2.set_title("Humanness Top-1 Match Rate", fontsize=10, fontweight="bold")
    ax2.legend(fontsize=8)

    fig.suptitle("Sapiens Humanness / Immunogenicity-Risk Summary", fontsize=11, fontweight="bold")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path


if __name__ == "__main__":
    main()
