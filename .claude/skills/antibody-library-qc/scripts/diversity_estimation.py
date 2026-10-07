#!/usr/bin/env python3
"""
Estimate the true total diversity (unique clone count) of a sequenced library
from a frequency-of-frequencies table, using two independent nonparametric
estimators: Chao1 and a Good-Turing-style Poisson coverage model.

This is a different question from qc_pipeline.py's "how clean is my library"
— it answers "how many *distinct* clones does my library actually contain,
including ones sequencing missed entirely?" It needs only a read-count column
(one row per observed unique clone/sequence, with how many reads it got) —
no design/PWM workbook required. Use it on the raw/unfiltered sequencing
output, not a design-rule-scored subset, since the estimators assume the
input is the full observed clone population.

Both estimators work off just two numbers beyond the raw totals: f1
(singletons — clones seen exactly once) and f2 (doubletons — seen exactly
twice). A library dominated by singletons (f1 >> f2) means sequencing barely
scratched the surface of true diversity, which is exactly when these
estimates diverge furthest from the naively-observed unique count.

Streams the input in chunks so multi-million-row / gzipped files don't need
to fit in memory at once.

Example:
    python diversity_estimation.py \\
        --input-csv clean_pairs_annotated.csv.gz \\
        --read-count-col read_count \\
        --output-dir diversity_output
"""

import argparse
import json
import math
import os
from collections import Counter

import pandas as pd

Z_SCORES = {
    "80%": 1.282,
    "90%": 1.645,
    "95%": 1.960,
    "99%": 2.576,
}


def stream_read_count_frequencies(file_path, read_count_col, chunksize=100_000):
    """Stream `file_path` (plain or gzip-compressed CSV, inferred from the
    extension) and return (freq_counter, total_unique_clones, total_reads),
    where freq_counter maps read-count value -> number of clones observed
    that many times."""
    freq_counter = Counter()
    total_unique_clones = 0
    total_reads = 0

    for chunk in pd.read_csv(file_path, usecols=[read_count_col], chunksize=chunksize):
        counts = chunk[read_count_col].dropna().astype(int)
        total_unique_clones += len(counts)
        total_reads += counts.sum()
        freq_counter.update(counts)

    return freq_counter, total_unique_clones, total_reads


def chao1_multi_ci(S_obs, f1, f2, z_scores=Z_SCORES):
    """Chao1 richness estimate plus log-scale confidence intervals at each
    level in `z_scores`. Falls back to the observed count (no extrapolation
    possible) when there are no doubletons or no singletons."""
    if f2 > 0 and f1 > 0:
        f1_f2_ratio = f1 / f2
        f_hat = (f1 ** 2) / (2.0 * f2)
        S_chao1 = S_obs + f_hat
        var_S = f2 * (((f1_f2_ratio / 4.0) ** 4) + (f1_f2_ratio ** 3) + ((f1_f2_ratio / 2.0) ** 2))
    else:
        return S_obs, {level: (S_obs, S_obs) for level in z_scores}

    cis = {}
    for level, z in z_scores.items():
        if var_S > 0 and f_hat > 0:
            C_factor = math.exp(z * math.sqrt(math.log(1.0 + (var_S / (f_hat ** 2)))))
            ci_lower = S_obs + (f_hat / C_factor)
            ci_upper = S_obs + (f_hat * C_factor)
        else:
            ci_lower, ci_upper = S_chao1, S_chao1
        cis[level] = (ci_lower, ci_upper)

    return S_chao1, cis


def poisson_solver(N, k):
    """Solve N/C = -ln(1 - k/C) for total diversity C, given N total reads
    and k unique clones observed, via bisection."""
    if N <= k or k <= 0:
        return float('inf')
    low, high = float(k), 1e15
    for _ in range(200):
        mid = (low + high) / 2.0
        if mid * (1.0 - math.exp(-N / mid)) < k:
            low = mid
        else:
            high = mid
    return (low + high) / 2.0


def poisson_multi_ci(N, k, z_scores=Z_SCORES):
    """Poisson coverage-model diversity estimate plus confidence intervals,
    propagated through poisson_solver from a normal approximation on the
    number of "duplicate" reads (N - k)."""
    S_poisson = poisson_solver(N, k)
    duplicates = N - k

    cis = {}
    if duplicates > 0:
        se_dup = math.sqrt(duplicates)
        for level, z in z_scores.items():
            k_upper = N - max(1.0, duplicates - z * se_dup)
            k_lower = N - (duplicates + z * se_dup)
            cis[level] = (poisson_solver(N, k_lower), poisson_solver(N, k_upper))
    else:
        cis = {level: (S_poisson, S_poisson) for level in z_scores}

    return S_poisson, cis


def estimate_diversity(file_path, read_count_col, chunksize=100_000):
    """End-to-end: stream the file, then return a dict with every number a
    caller would want (point estimates, CIs, and the raw f1/f2/N/S_obs
    inputs they were derived from)."""
    freq_counter, S_obs, N = stream_read_count_frequencies(file_path, read_count_col, chunksize)
    f1 = freq_counter.get(1, 0)
    f2 = freq_counter.get(2, 0)

    s_chao, chao_cis = chao1_multi_ci(float(S_obs), f1, f2)
    s_pois, pois_cis = poisson_multi_ci(float(N), float(S_obs))

    return {
        "total_reads": int(N),
        "unique_clones_observed": int(S_obs),
        "singletons_f1": int(f1),
        "doubletons_f2": int(f2),
        "chao1_estimate": s_chao,
        "poisson_estimate": s_pois,
        "confidence_intervals": {
            level: {
                "chao1_lower": chao_cis[level][0],
                "chao1_upper": chao_cis[level][1],
                "poisson_lower": pois_cis[level][0],
                "poisson_upper": pois_cis[level][1],
            }
            for level in Z_SCORES
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input-csv", required=True,
                         help="CSV or CSV.GZ with one row per unique observed clone/sequence.")
    parser.add_argument("--read-count-col", default="read_count",
                         help="Column giving each unique clone's read count.")
    parser.add_argument("--chunksize", type=int, default=100_000)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print(f"Streaming '{args.read_count_col}' from {args.input_csv} in chunks of {args.chunksize:,}...")
    result = estimate_diversity(args.input_csv, args.read_count_col, args.chunksize)

    print("\n=== DATASET SUMMARY ===")
    print(f"Total Reads (N)         : {result['total_reads']:,}")
    print(f"Unique Clones (S_obs)   : {result['unique_clones_observed']:,}")
    print(f"Singletons (f1, count=1): {result['singletons_f1']:,}")
    print(f"Doubletons (f2, count=2): {result['doubletons_f2']:,}")

    print("\n=== POINT ESTIMATES ===")
    print(f"Chao1 Estimate  : {result['chao1_estimate']:,.0f} unique clones")
    print(f"Poisson Estimate: {result['poisson_estimate']:,.0f} unique clones")

    print("\n=== CONFIDENCE INTERVAL TABLE ===")
    rows = []
    for level, ci in result["confidence_intervals"].items():
        rows.append({
            "Confidence Level": level,
            "Chao1 Lower CI": f"{ci['chao1_lower']:,.0f}",
            "Chao1 Upper CI": f"{ci['chao1_upper']:,.0f}",
            "Poisson Lower CI": f"{ci['poisson_lower']:,.0f}",
            "Poisson Upper CI": f"{ci['poisson_upper']:,.0f}",
        })
    print(pd.DataFrame(rows).to_string(index=False))

    output_path = os.path.join(args.output_dir, "diversity_estimate.json")
    with open(output_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nSaved {output_path}")


if __name__ == "__main__":
    main()
