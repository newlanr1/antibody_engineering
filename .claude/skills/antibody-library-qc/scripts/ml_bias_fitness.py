#!/usr/bin/env python3
"""
Run after qc_pipeline.py. Computes KL-divergence design bias per designed
position (always, on the full dataset — it's cheap) and protein-LM pseudo-
perplexity fitness scoring (on by default too, but automatically capped to a
random representative sample for large files, since it's a transformer
forward pass per sequence and much slower per-read than the alignment work
in qc_pipeline.py).

Fitness scoring defaults to AbLang2, an antibody-specific LM (vs. the earlier
default, ESM-2, trained on all of UniRef — general proteins, not specifically
antibodies). AbLang2 scores heavy+light jointly; this script uses it to score
each chain alone (preserving the existing per-chain plots) *and* adds a joint
"pairing compatibility" score as an extra metric when both chains are given.
Pass --fitness-backend esm2 to use the old general-purpose model instead.

Example:
    python ml_bias_fitness.py \\
        --scored-reads-csv antibody_qc_output/scored_reads.csv \\
        --design-xlsx AffMat_design.xlsx \\
        --heavy-sheet "Heavy Chain Design" --heavy-col HAA \\
        --light-sheet "Light Chain Design" --light-col LAA \\
        --read-count-col read_count \\
        --output-dir antibody_qc_output

Add --skip-fitness-lm to run KL divergence only, or --fitness-sample-size 0
to force scoring every row regardless of size (expect this to take hours
past roughly 100k rows on CPU, for either backend).
"""

import argparse
import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pwm_extraction import extract_full_pwm
from bias_and_fitness import (
    calculate_pwm_kl_divergence,
    score_esm2_perplexity_batch,
    compute_wt_reference_perplexity,
    load_ablang2,
    score_ablang2_chain_batch,
    score_ablang2_paired_batch,
    compute_ablang2_wt_reference,
    compute_ablang2_paired_wt_reference,
    plot_perplexity_vs_errors,
)

# AbLang2's single-forward-pass scoring call runs ~5-6x slower per sequence
# than ESM-2's (benchmarked ~45ms vs ~7ms), and when both chains are present
# this script makes 3 calls per row (heavy-alone, light-alone, paired)
# instead of ESM-2's 2 (one per chain). A smaller default sample keeps
# end-to-end runtime in the same ballpark as the old ESM-2 default.
DEFAULT_SAMPLE_SIZE = {"ablang2": 15_000, "esm2": 50_000}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--scored-reads-csv", required=True,
                         help="scored_reads.csv from qc_pipeline.py (must contain the "
                              "{Label}_Design/{Label}_Intended_Mutations/{Label}_True_Errors "
                              "columns it writes).")
    parser.add_argument("--design-xlsx", required=True)
    parser.add_argument("--read-count-col", default="read_count")
    parser.add_argument("--heavy-sheet", default=None)
    parser.add_argument("--heavy-col", default=None)
    parser.add_argument("--light-sheet", default=None)
    parser.add_argument("--light-col", default=None)
    parser.add_argument("--error-threshold", type=int, default=2,
                         help="Same meaning as in qc_pipeline.py — used only to draw the "
                              "QC cutoff line on the perplexity scatter plot.")
    parser.add_argument("--fitness-backend", choices=["ablang2", "esm2"], default="ablang2",
                         help="Which protein LM to use for fitness scoring. ablang2 (default) "
                              "is antibody-specific; esm2 is the earlier general-purpose model, "
                              "kept available for comparison.")
    parser.add_argument("--esm2-model", default="facebook/esm2_t6_8M_UR50D",
                         help="Hugging Face model id, only used when --fitness-backend esm2.")
    parser.add_argument("--skip-fitness-lm", action="store_true",
                         help="Skip fitness-LM scoring entirely (KL divergence only).")
    parser.add_argument("--fitness-sample-size", type=int, default=None,
                         help="If the dataset has more rows than this, fitness scoring runs on "
                              "a random sample of this many rows instead of everything. Defaults "
                              f"to {DEFAULT_SAMPLE_SIZE['ablang2']:,} for ablang2 or "
                              f"{DEFAULT_SAMPLE_SIZE['esm2']:,} for esm2, reflecting their "
                              "different per-sequence costs. Pass 0 to force scoring every row.")
    parser.add_argument("--fitness-sample-seed", type=int, default=42,
                         help="Random seed for the fitness-LM sample, so repeat runs are reproducible.")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    if not args.heavy_sheet and not args.light_sheet:
        parser.error("Provide at least --heavy-sheet/--heavy-col.")

    sample_size = args.fitness_sample_size
    if sample_size is None:
        sample_size = DEFAULT_SAMPLE_SIZE[args.fitness_backend]

    os.makedirs(args.output_dir, exist_ok=True)
    df = pd.read_csv(args.scored_reads_csv)

    chains = []
    if args.heavy_sheet:
        rules = extract_full_pwm(args.design_xlsx, args.heavy_sheet)
        chains.append(dict(label="Heavy", side="heavy", seq_col=args.heavy_col, design_col="Heavy_Design",
                            intended_col="Heavy_Intended_Mutations", error_col="Heavy_True_Errors",
                            window_col="Heavy_Best_Window", rules=rules))
    if args.light_sheet:
        rules = extract_full_pwm(args.design_xlsx, args.light_sheet)
        chains.append(dict(label="Light", side="light", seq_col=args.light_col, design_col="Light_Design",
                            intended_col="Light_Intended_Mutations", error_col="Light_True_Errors",
                            window_col="Light_Best_Window", rules=rules))

    # --- KL divergence (always, it's cheap — especially now that it reuses the
    # alignment windows qc_pipeline.py already computed, instead of re-aligning) ---
    kl_frames = []
    for chain in chains:
        print(f"Calculating KL divergence for {chain['label']} chain designs...")
        window_col = chain["window_col"] if chain["window_col"] in df.columns else None
        kl_frames.append(calculate_pwm_kl_divergence(
            df, chain["rules"], chain["seq_col"], chain["design_col"],
            args.read_count_col, chain["label"], best_window_col=window_col,
        ))
    kl_df = pd.concat(kl_frames, ignore_index=True) if kl_frames else pd.DataFrame()
    kl_csv_path = os.path.join(args.output_dir, "kl_divergence.csv")
    kl_df.to_csv(kl_csv_path, index=False)
    print(f"Saved {kl_csv_path}")

    if not kl_df.empty:
        print("\nTop 10 positions by design-bias (highest KL divergence):")
        print(kl_df.sort_values(by='KL_Divergence', ascending=False).head(10).to_string(index=False))

    # --- Fitness LM scoring (on by default, auto-sampled for large files) ---
    perplexity_images = []
    fitness_sampled = False
    fitness_rows_scored = 0
    model_label = "AbLang2" if args.fitness_backend == "ablang2" else "ESM-2"

    if args.skip_fitness_lm:
        print("--skip-fitness-lm given — skipping fitness-LM scoring.")
    else:
        fitness_df = df
        if sample_size and len(df) > sample_size:
            fitness_df = df.sample(n=sample_size, random_state=args.fitness_sample_seed)
            fitness_sampled = True
        fitness_rows_scored = len(fitness_df)
        fitness_df = fitness_df.copy()
        sample_note = f" (random sample of {fitness_rows_scored:,} of {len(df):,} rows)" if fitness_sampled else ""

        # Score the *aligned* window (qc_pipeline.py's {Label}_Best_Window), not the raw
        # sequence column. Raw reads routinely carry CDR3/J-segment/constant-domain
        # carryover that the design's wt_sequence (used for the WT fold-compatible
        # reference below) never includes -- scoring the untrimmed read against a
        # baseline built from the trimmed WT would inflate every read's perplexity
        # relative to WT regardless of actual mutation/error count, which is exactly
        # the kind of apples-to-oranges comparison this metric is supposed to avoid.
        # Falls back to the raw column only if a window wasn't available (standalone
        # use without qc_pipeline.py having run first).
        for chain in chains:
            chain["score_col"] = chain["window_col"] if chain["window_col"] in fitness_df.columns else chain["seq_col"]

        if args.fitness_backend == "ablang2":
            print(f"Loading AbLang2...")
            ablang_model = load_ablang2()

            for chain in chains:
                print(f"Scoring {chain['label']} chain sequences with AbLang2{sample_note}...")
                scores = score_ablang2_chain_batch(fitness_df[chain["score_col"]], ablang_model, chain["side"])
                perplexity_col = f"{chain['label']}_Fitness_Perplexity"
                fitness_df[perplexity_col] = fitness_df[chain["score_col"]].map(scores)

                fold_reference = compute_ablang2_wt_reference(chain["rules"], ablang_model, chain["side"])

                plot_path = os.path.join(args.output_dir, f"{chain['label'].lower()}_fitness_perplexity_vs_errors.png")
                plot_perplexity_vs_errors(fitness_df, chain["error_col"], perplexity_col, chain["intended_col"],
                                           f"{chain['label']} Chain", plot_path, error_threshold=args.error_threshold,
                                           fold_reference=fold_reference, model_label=model_label)
                perplexity_images.append([f"{chain['label']} Chain: AbLang2 Fitness vs. Framework Errors", plot_path])
                print(f"Saved {plot_path} (WT baseline: mean {fold_reference['mean']:.3f}, range {fold_reference['min']:.3f}-{fold_reference['max']:.3f})")

            # Joint pairing-compatibility score, only meaningful with both chains.
            if len(chains) == 2:
                heavy, light = chains
                print(f"Scoring Heavy+Light pairing compatibility with AbLang2{sample_note}...")
                pair_scores = score_ablang2_paired_batch(
                    fitness_df[heavy["score_col"]], fitness_df[light["score_col"]], ablang_model)
                fitness_df["Paired_Fitness_Perplexity"] = [
                    pair_scores.get((h, l), float("nan"))
                    for h, l in zip(fitness_df[heavy["score_col"]], fitness_df[light["score_col"]])
                ]
                fitness_df["Combined_True_Errors"] = fitness_df[heavy["error_col"]] + fitness_df[light["error_col"]]
                fitness_df["Combined_Intended_Mutations"] = (
                    fitness_df[heavy["intended_col"]] + fitness_df[light["intended_col"]])

                paired_fold_reference = compute_ablang2_paired_wt_reference(
                    heavy["rules"], light["rules"], ablang_model)

                plot_path = os.path.join(args.output_dir, "paired_fitness_perplexity_vs_errors.png")
                plot_perplexity_vs_errors(fitness_df, "Combined_True_Errors", "Paired_Fitness_Perplexity",
                                           "Combined_Intended_Mutations", "Heavy+Light Pair", plot_path,
                                           error_threshold=args.error_threshold * 2,
                                           fold_reference=paired_fold_reference, model_label=model_label)
                perplexity_images.append(["Heavy+Light Pair: AbLang2 Pairing Compatibility vs. Framework Errors",
                                           plot_path])
                print(f"Saved {plot_path} (WT-pair baseline: mean {paired_fold_reference['mean']:.3f}, range {paired_fold_reference['min']:.3f}-{paired_fold_reference['max']:.3f})")

        else:  # esm2
            for chain in chains:
                print(f"Scoring {chain['label']} chain sequences with ESM-2 ({args.esm2_model}){sample_note}...")
                scores = score_esm2_perplexity_batch(fitness_df[chain["score_col"]], model_name=args.esm2_model)
                perplexity_col = f"{chain['label']}_Fitness_Perplexity"
                fitness_df[perplexity_col] = fitness_df[chain["score_col"]].map(scores)

                fold_reference = compute_wt_reference_perplexity(chain["rules"], model_name=args.esm2_model)

                plot_path = os.path.join(args.output_dir, f"{chain['label'].lower()}_fitness_perplexity_vs_errors.png")
                plot_perplexity_vs_errors(fitness_df, chain["error_col"], perplexity_col, chain["intended_col"],
                                           f"{chain['label']} Chain", plot_path, error_threshold=args.error_threshold,
                                           fold_reference=fold_reference, model_label=model_label)
                perplexity_images.append([f"{chain['label']} Chain: ESM-2 Fitness vs. Framework Errors", plot_path])
                print(f"Saved {plot_path} (WT baseline: mean {fold_reference['mean']:.3f}, range {fold_reference['min']:.3f}-{fold_reference['max']:.3f})")

        scored_with_fitness_path = os.path.join(
            args.output_dir, "scored_reads_with_fitness_sample.csv" if fitness_sampled else "scored_reads_with_fitness.csv")
        fitness_df.to_csv(scored_with_fitness_path, index=False)
        print(f"Saved {scored_with_fitness_path}")

    manifest_path = os.path.join(args.output_dir, "ml_manifest.json")
    with open(manifest_path, "w") as f:
        json.dump({
            "kl_divergence_csv": kl_csv_path,
            "perplexity_images": perplexity_images,
            "fitness_backend": args.fitness_backend,
            "fitness_sampled": fitness_sampled,
            "fitness_rows_scored": fitness_rows_scored,
            "fitness_total_rows": len(df),
        }, f, indent=2)
    print(f"ML manifest: {manifest_path}")


if __name__ == "__main__":
    main()
