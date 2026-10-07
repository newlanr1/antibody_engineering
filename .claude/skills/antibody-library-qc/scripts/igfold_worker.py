#!/usr/bin/env python3
"""
Standalone worker: predicts a 3D structure for each paired heavy+light
sequence in --input-csv using IgFold. Must be run with the DEDICATED IgFold
venv's interpreter (Python >=3.11) -- this is a different Python version
than the rest of this skill, so it can't be imported directly; it's invoked
via subprocess from structure_prediction.py in the main venv.

IgFold's code/weights are licensed for non-commercial use only (JHU Academic
Software License) -- commercial use needs a separate license through JHU
Technology Ventures. See https://github.com/Graylab/IgFold/blob/main/LICENSE.md.

Runs with do_refine=False (no OpenMM/energy-minimization step) for speed and
to avoid an extra heavy dependency -- benchmarked at ~1.2s/structure on CPU
with a single model, vs. minutes+ for a refined multi-model ensemble. That's
fast enough that this worker could in principle run on more than a handful
of sequences, but structure_prediction.py still only ever calls it on a
small, explicit shortlist (see that script) -- this is still 1000x+ slower
per sequence than the alignment work the rest of the skill does, and nothing
here has been checked against a large-scale run.

Example:
    python igfold_worker.py \\
        --input-csv shortlist.csv --id-col id --heavy-col heavy_seq --light-col light_seq \\
        --output-dir structures/ --results-csv structure_results.csv
"""

import argparse
import os

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input-csv", required=True)
    parser.add_argument("--id-col", default="id")
    parser.add_argument("--heavy-col", default="heavy_seq")
    parser.add_argument("--light-col", default="light_seq")
    parser.add_argument("--output-dir", required=True, help="Directory to write one PDB per row into.")
    parser.add_argument("--results-csv", required=True, help="Where to write id/pdb_path/mean_prmsd results.")
    args = parser.parse_args()

    from igfold import IgFoldRunner

    os.makedirs(args.output_dir, exist_ok=True)
    df = pd.read_csv(args.input_csv)

    print(f"Loading IgFold ({len(df)} structures to predict)...")
    igfold = IgFoldRunner(num_models=1, try_gpu=False)

    results = []
    for _, row in df.iterrows():
        seq_id = str(row[args.id_col])
        heavy = row[args.heavy_col]
        light = row[args.light_col]
        pdb_path = os.path.join(args.output_dir, f"{seq_id}.pdb")

        out = igfold.fold(pdb_path, sequences={"H": heavy, "L": light}, do_refine=False, do_renum=False)
        mean_prmsd = float(out.prmsd.mean().item())

        results.append({"id": seq_id, "pdb_path": pdb_path, "mean_predicted_rmsd": round(mean_prmsd, 3)})
        print(f"  {seq_id}: {pdb_path} (mean predicted RMSD {mean_prmsd:.3f} Å)")

    pd.DataFrame(results).to_csv(args.results_csv, index=False)
    print(f"Saved {args.results_csv}")


if __name__ == "__main__":
    main()
