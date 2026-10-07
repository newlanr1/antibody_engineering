#!/usr/bin/env python3
"""
Standalone worker: scores each sequence in --input-csv with Sapiens (a human
antibody language model) for humanness/immunogenicity-risk signal. Must be
run with the DEDICATED Sapiens environment's interpreter (Python 3.7/3.8 --
newer Python breaks Sapiens per its own README) -- invoked via subprocess
from humanness_scoring.py in the main venv.

For each sequence, Sapiens' predict_scores gives a per-position probability
over the 20 amino acids (what a human antibody would have at that position,
per the model). Two metrics are derived from that: a perplexity-style score
(exp(mean(-log(prob of the observed residue))), consistent framing with the
fitness-LM metrics elsewhere in this skill -- lower = more human-like) and a
simpler "% positions matching the model's top-1 human residue", closer to
how OASis/Sapiens humanness is usually reported in the literature.

Example:
    python sapiens_worker.py \\
        --input-csv sample.csv --seq-col sequence --chain-col chain --output-csv results.csv
"""

import argparse
import math

import numpy as np
import pandas as pd
import sapiens


def score_sequence(sequence, chain):
    if not isinstance(sequence, str) or not sequence:
        return None, None
    scores = sapiens.predict_scores(sequence, chain)  # (len(sequence), 20) amino-acid probabilities
    observed_probs = []
    match_count = 0
    for pos, residue in enumerate(sequence):
        if pos >= len(scores) or residue not in scores.columns:
            continue
        row = scores.iloc[pos]
        prob = row.get(residue, 1e-6)
        observed_probs.append(max(prob, 1e-6))
        if row.idxmax() == residue:
            match_count += 1
    if not observed_probs:
        return None, None
    perplexity = round(math.exp(-np.mean(np.log(observed_probs))), 3)
    match_pct = round(100 * match_count / len(observed_probs), 2)
    return perplexity, match_pct


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input-csv", required=True,
                         help="CSV with a sequence column and a chain column (values 'H' or 'L').")
    parser.add_argument("--seq-col", default="sequence")
    parser.add_argument("--chain-col", default="chain")
    parser.add_argument("--output-csv", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.input_csv)
    print(f"Scoring {len(df)} sequences with Sapiens...")

    perplexities, match_pcts = [], []
    for seq, chain in zip(df[args.seq_col], df[args.chain_col]):
        ppl, match_pct = score_sequence(seq, chain)
        perplexities.append(ppl)
        match_pcts.append(match_pct)

    df["humanness_perplexity"] = perplexities
    df["humanness_match_pct"] = match_pcts
    df.to_csv(args.output_csv, index=False)
    print(f"Saved {args.output_csv}")


if __name__ == "__main__":
    main()
