#!/usr/bin/env python3
"""
Score an antibody library NGS dataset against a design-rule workbook and
produce QC metrics + plots. Run this first; it does not write a PowerPoint —
see build_report.py for that, which needs real conclusions as input.

Supports one chain (e.g. a single-domain/nanobody/scFv library) or two chains
(paired Heavy + Light). Chain column names, sheet names, and thresholds are
all parameters — nothing here is hardcoded to a specific dataset.

Example (paired heavy/light library):
    python qc_pipeline.py \\
        --sequences-csv preview_1000_pairs.csv \\
        --design-xlsx AffMat_design.xlsx \\
        --heavy-sheet "Heavy Chain Design" --heavy-col HAA \\
        --light-sheet "Light Chain Design" --light-col LAA \\
        --read-count-col read_count \\
        --output-dir antibody_qc_output

Example (single chain library):
    python qc_pipeline.py \\
        --sequences-csv nanobody_reads.csv \\
        --design-xlsx Nb_design.xlsx \\
        --heavy-sheet "Nanobody Design" --heavy-col AA \\
        --read-count-col read_count \\
        --output-dir nanobody_qc_output
"""

import argparse
import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pwm_extraction import extract_full_pwm
from aligner import score_reads_parallel
from plotting import plot_error_distributions, generate_split_observed_pwm


def score_chain(df, chain_label, seq_col, rules_dict, n_processes=None):
    """Score every read in df[seq_col] against rules_dict (in parallel — this
    is the single most expensive step in the whole pipeline, since it tests
    each read against every design) and write the results, including the
    winning alignment window, back onto df. Returns the new column names so
    callers don't have to hardcode the f"{chain_label}_..." pattern."""
    design_col = f"{chain_label}_Design"
    intended_col = f"{chain_label}_Intended_Mutations"
    error_col = f"{chain_label}_True_Errors"
    window_col = f"{chain_label}_Best_Window"

    results = score_reads_parallel(df[seq_col], rules_dict, n_processes=n_processes)
    # NB: each result is a pd.Series with a default int index (0..3), so
    # pd.DataFrame(results, columns=[...]) would try to *select* those string
    # labels from that int index and silently come back all-NaN. Build with
    # the default int columns first, then rename positionally.
    results_df = pd.DataFrame(results)
    results_df.columns = [design_col, intended_col, error_col, window_col]
    results_df.index = df.index
    df[[design_col, intended_col, error_col, window_col]] = results_df
    return design_col, intended_col, error_col, window_col


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sequences-csv", required=True,
                         help="CSV of observed reads with amino acid sequence column(s) "
                              "and a read-count/abundance column.")
    parser.add_argument("--design-xlsx", required=True,
                         help="Excel workbook with per-gene wild-type sequence + allowed "
                              "mutation design rules (see pwm_extraction.py for expected layout).")
    parser.add_argument("--read-count-col", default="read_count",
                         help="Column in sequences-csv giving each row's read abundance.")
    parser.add_argument("--heavy-sheet", default=None, help="Sheet name for the first/heavy chain.")
    parser.add_argument("--heavy-col", default=None, help="Sequence column name for the first/heavy chain.")
    parser.add_argument("--light-sheet", default=None, help="Sheet name for the second/light chain.")
    parser.add_argument("--light-col", default=None, help="Sequence column name for the second/light chain.")
    parser.add_argument("--error-threshold", type=int, default=2,
                         help="Max true framework errors for a read to count as 'valid'.")
    parser.add_argument("--fidelity-threshold", type=float, default=0.85,
                         help="Minimum observed wild-type frequency at a framework position "
                              "before it's flagged as a hotspot in the PWM plots.")
    parser.add_argument("--output-dir", required=True, help="Directory to write outputs into.")
    parser.add_argument("--n-processes", type=int, default=None,
                         help="CPU cores to parallelize read scoring across. Defaults to all "
                              "available cores (os.cpu_count()). Scoring is embarrassingly "
                              "parallel per-read, so this is the main lever for large files — "
                              "set to 1 to force single-process (e.g. for debugging).")
    args = parser.parse_args()

    if not args.heavy_sheet and not args.light_sheet:
        parser.error("Provide at least --heavy-sheet/--heavy-col (single-chain libraries use only this pair).")
    if bool(args.heavy_sheet) != bool(args.heavy_col):
        parser.error("--heavy-sheet and --heavy-col must be given together.")
    if bool(args.light_sheet) != bool(args.light_col):
        parser.error("--light-sheet and --light-col must be given together.")

    os.makedirs(args.output_dir, exist_ok=True)

    # Only load the columns this pipeline actually touches. Real NGS exports often
    # carry dozens of annotation columns (allele calls, lengths, provenance, ...)
    # that would otherwise multiply memory use and scored_reads.csv's size for no
    # benefit — irrelevant at 1,000 rows, but matters once you're at millions.
    usecols = [c for c in [args.heavy_col, args.light_col, args.read_count_col] if c]
    df = pd.read_csv(args.sequences_csv, usecols=usecols)
    total_reads = df[args.read_count_col].sum()
    total_rows = len(df)

    chains = []  # list of dicts: label, seq_col, design_col, intended_col, error_col, window_col, rules

    if args.heavy_sheet:
        print(f"Loading design rules from '{args.heavy_sheet}'...")
        heavy_rules = extract_full_pwm(args.design_xlsx, args.heavy_sheet)
        print(f"  {len(heavy_rules)} design(s) loaded.")
        print(f"Scoring {total_rows:,} reads against heavy-chain designs "
              f"(across {args.n_processes or os.cpu_count() or 1} process(es))...")
        design_col, intended_col, error_col, window_col = score_chain(
            df, "Heavy", args.heavy_col, heavy_rules, n_processes=args.n_processes)
        chains.append(dict(label="Heavy", seq_col=args.heavy_col, design_col=design_col,
                            intended_col=intended_col, error_col=error_col,
                            window_col=window_col, rules=heavy_rules))

    if args.light_sheet:
        print(f"Loading design rules from '{args.light_sheet}'...")
        light_rules = extract_full_pwm(args.design_xlsx, args.light_sheet)
        print(f"  {len(light_rules)} design(s) loaded.")
        print(f"Scoring {total_rows:,} reads against light-chain designs "
              f"(across {args.n_processes or os.cpu_count() or 1} process(es))...")
        design_col, intended_col, error_col, window_col = score_chain(
            df, "Light", args.light_col, light_rules, n_processes=args.n_processes)
        chains.append(dict(label="Light", seq_col=args.light_col, design_col=design_col,
                            intended_col=intended_col, error_col=error_col,
                            window_col=window_col, rules=light_rules))

    # --- Metrics ---
    summary_metrics = {
        "Total Sequence Rows Analyzed": f"{total_rows:,}",
        "Total Read Abundance Sum": f"{total_reads:,}",
    }

    valid_masks = []
    for chain in chains:
        valid_mask = df[chain["error_col"]] <= args.error_threshold
        valid_masks.append(valid_mask)
        valid_reads = df.loc[valid_mask, args.read_count_col].sum()
        efficiency = (valid_reads / total_reads) * 100 if total_reads else 0.0
        summary_metrics[f"{chain['label']} Chain Valid Reads (≤{args.error_threshold} Errors)"] = \
            f"{valid_reads:,} ({efficiency:.1f}%)"

    if len(chains) == 2:
        both_valid = valid_masks[0] & valid_masks[1]
        valid_both = df[both_valid]
        paired_efficiency = (valid_both[args.read_count_col].sum() / total_reads) * 100 if total_reads else 0.0
        summary_metrics["Fully Functional Paired Reads"] = \
            f"{valid_both[args.read_count_col].sum():,} ({paired_efficiency:.1f}%)"

        seq_cols = [c["seq_col"] for c in chains]
        unique_raw_pairs = df.groupby(seq_cols)[args.read_count_col].sum().reset_index()
        unique_valid_pairs = valid_both.groupby(seq_cols)[args.read_count_col].sum().reset_index()
        summary_metrics["Total Unique Sequence Pairs"] = f"{len(unique_raw_pairs):,}"
        summary_metrics["Unique Quality-Filtered Pairs"] = f"{len(unique_valid_pairs):,}"

    print("\n--- LIBRARY QC SUMMARY METRICS ---")
    for k, v in summary_metrics.items():
        print(f"{k}: {v}")

    # --- Plots ---
    distribution_images = []
    pwm_images = []

    for chain in chains:
        dist_path = os.path.join(args.output_dir, f"{chain['label'].lower()}_chain_distributions.png")
        plot_error_distributions(df, chain["intended_col"], chain["error_col"], args.read_count_col,
                                  f"{chain['label']} Chain", dist_path)
        distribution_images.append([f"{chain['label']} Chain Mutagenesis & Error Profiles", dist_path])
        print(f"Saved {dist_path}")

        pwm_dir = os.path.join(args.output_dir, "split_pwm_plots")
        paths_by_gene = generate_split_observed_pwm(
            df, chain["rules"], chain["seq_col"], chain["design_col"], args.read_count_col,
            chain["label"], pwm_dir, fidelity_threshold=args.fidelity_threshold,
            best_window_col=chain["window_col"],
        )
        for gene, path in paths_by_gene.items():
            pwm_images.append([f"{chain['label']} Chain PWM (FW vs CDR): {gene}", path])
        print(f"Saved {len(paths_by_gene)} split PWM plot(s) for {chain['label']} chain to {pwm_dir}/")

    # --- Persist machine-readable outputs for the next step (build_report.py) ---
    scored_csv_path = os.path.join(args.output_dir, "scored_reads.csv")
    df.to_csv(scored_csv_path, index=False)

    metrics_path = os.path.join(args.output_dir, "qc_metrics.json")
    with open(metrics_path, "w") as f:
        json.dump(summary_metrics, f, indent=2)

    manifest_path = os.path.join(args.output_dir, "plot_manifest.json")
    with open(manifest_path, "w") as f:
        json.dump({"distribution_images": distribution_images, "pwm_images": pwm_images}, f, indent=2)

    print(f"\nScored reads: {scored_csv_path}")
    print(f"Metrics JSON: {metrics_path}")
    print(f"Plot manifest: {manifest_path}")


if __name__ == "__main__":
    main()
