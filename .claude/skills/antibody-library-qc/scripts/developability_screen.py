#!/usr/bin/env python3
"""
Run after qc_pipeline.py. Screens every read for common sequence-based
developability liabilities (deamidation, isomerization, oxidation-prone CDR
residues, N-glycosylation motifs, cysteine-count anomalies, CDR hydrophobic
patches). See developability.py's module docstring for why this is a
heuristic sequence screen rather than a trained ML model, and why it runs on
the full dataset (cheap -- no model, no re-alignment) rather than a sample.

Example:
    python developability_screen.py \\
        --scored-reads-csv antibody_qc_output/scored_reads.csv \\
        --design-xlsx AffMat_design.xlsx \\
        --heavy-sheet "Heavy Chain Design" --heavy-col HAA \\
        --light-sheet "Light Chain Design" --light-col LAA \\
        --read-count-col read_count \\
        --output-dir antibody_qc_output
"""

import argparse
import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pwm_extraction import extract_full_pwm
from developability import (
    score_developability_batch,
    summarize_liability_prevalence,
    plot_liability_prevalence,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--scored-reads-csv", required=True,
                         help="scored_reads.csv from qc_pipeline.py (must contain the "
                              "{Label}_Design/{Label}_Best_Window columns it writes).")
    parser.add_argument("--design-xlsx", required=True)
    parser.add_argument("--read-count-col", default="read_count")
    parser.add_argument("--heavy-sheet", default=None)
    parser.add_argument("--heavy-col", default=None)
    parser.add_argument("--light-sheet", default=None)
    parser.add_argument("--light-col", default=None)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    if not args.heavy_sheet and not args.light_sheet:
        parser.error("Provide at least --heavy-sheet/--heavy-col.")

    os.makedirs(args.output_dir, exist_ok=True)
    df = pd.read_csv(args.scored_reads_csv)

    chains = []
    if args.heavy_sheet:
        rules = extract_full_pwm(args.design_xlsx, args.heavy_sheet)
        chains.append(dict(label="Heavy", design_col="Heavy_Design", window_col="Heavy_Best_Window", rules=rules))
    if args.light_sheet:
        rules = extract_full_pwm(args.design_xlsx, args.light_sheet)
        chains.append(dict(label="Light", design_col="Light_Design", window_col="Light_Best_Window", rules=rules))

    all_flags = []
    prevalence_by_chain = {}
    for chain in chains:
        print(f"Screening {chain['label']} chain reads for developability liabilities...")
        flags_df = score_developability_batch(df, chain["design_col"], chain["window_col"],
                                               chain["rules"], chain["label"])
        all_flags.append(flags_df)
        prevalence_by_chain[chain["label"]] = summarize_liability_prevalence(
            flags_df, df[args.read_count_col], chain["label"])

    flags_df = pd.concat(all_flags, axis=1) if all_flags else pd.DataFrame(index=df.index)
    flags_csv_path = os.path.join(args.output_dir, "developability_flags.csv")
    flags_df.to_csv(flags_csv_path, index=False)
    print(f"Saved {flags_csv_path}")

    print("\n--- DEVELOPABILITY LIABILITY PREVALENCE (abundance-weighted %) ---")
    for chain_label, prevalence in prevalence_by_chain.items():
        print(f"{chain_label}:")
        for label, pct in prevalence.items():
            print(f"  {label}: {pct}%")

    summary_path = os.path.join(args.output_dir, "developability_summary.json")
    with open(summary_path, "w") as f:
        json.dump(prevalence_by_chain, f, indent=2)
    print(f"\nSaved {summary_path}")

    chart_path = os.path.join(args.output_dir, "developability_liability_prevalence.png")
    plot_liability_prevalence(prevalence_by_chain, chart_path)
    print(f"Saved {chart_path}")

    manifest_path = os.path.join(args.output_dir, "developability_manifest.json")
    with open(manifest_path, "w") as f:
        json.dump({
            "flags_csv": flags_csv_path,
            "summary_json": summary_path,
            "prevalence_image": chart_path,
        }, f, indent=2)
    print(f"Developability manifest: {manifest_path}")


if __name__ == "__main__":
    main()
