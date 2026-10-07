#!/usr/bin/env python3
"""
Assemble the final QC PowerPoint from qc_pipeline.py's outputs plus
conclusions written after actually looking at the metrics and plots.

Methodology and conclusion bullets are passed in as JSON files (one JSON
array of strings each) rather than typed on the command line, since they're
usually multi-sentence findings. Write them with the Write tool, then run:

    python build_report.py \\
        --output-dir antibody_qc_output \\
        --title "Antibody Library Analysis" \\
        --subtitle "Framework & CDR Split PWM Profiling" \\
        --methodology-json methodology.json \\
        --conclusions-json conclusions.json \\
        --diversity-json diversity_output/diversity_estimate.json \\
        --pptx-name Library_QC_Summary.pptx
"""

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pptx_report import build_qc_deck


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-dir", required=True,
                         help="Directory containing qc_metrics.json and plot_manifest.json "
                              "from qc_pipeline.py. The .pptx is also written here.")
    parser.add_argument("--title", default="Antibody Library Analysis")
    parser.add_argument("--subtitle", default="NGS Library QC Summary")
    parser.add_argument("--methodology-json", required=True,
                         help="Path to a JSON array of methodology bullet strings.")
    parser.add_argument("--conclusions-json", default=None,
                         help="Path to a JSON array of conclusion bullet strings, based on "
                              "what this run's metrics/plots actually show. Omit if you don't "
                              "have real findings to report yet.")
    parser.add_argument("--ml-manifest", default=None,
                         help="Path to ml_manifest.json from ml_bias_fitness.py, if that step "
                              "was run. Adds its perplexity scatter plots as slides.")
    parser.add_argument("--ml-bullets-json", default=None,
                         help="Path to a JSON array of bullet strings summarizing the KL-divergence "
                              "/ ESM-2 findings actually observed in this run. Only used if "
                              "--ml-manifest is also given.")
    parser.add_argument("--diversity-json", default=None,
                         help="Path to diversity_estimate.json from diversity_estimation.py. "
                              "Adds a diversity summary slide + confidence-interval table slide. "
                              "This is a separate analysis from the rest of the deck (no design "
                              "workbook needed), so pass whatever path it was written to — it "
                              "doesn't need to live under --output-dir.")
    parser.add_argument("--developability-manifest", default=None,
                         help="Path to developability_manifest.json from developability_screen.py, "
                              "if that step was run. Adds the liability-prevalence bar chart as a slide.")
    parser.add_argument("--developability-bullets-json", default=None,
                         help="Path to a JSON array of bullet strings summarizing the liability-screen "
                              "findings actually observed in this run. Only used if "
                              "--developability-manifest is also given.")
    parser.add_argument("--humanness-manifest", default=None,
                         help="Path to humanness_manifest.json from humanness_scoring.py, if that "
                              "step was run. Adds the humanness summary bar chart as a slide. "
                              "No-ops if the manifest says skipped.")
    parser.add_argument("--humanness-bullets-json", default=None,
                         help="Path to a JSON array of bullet strings summarizing the humanness/ "
                              "immunogenicity findings actually observed in this run. Only used "
                              "if --humanness-manifest is also given.")
    parser.add_argument("--structure-manifest", default=None,
                         help="Path to structure_manifest.json from structure_prediction.py, if "
                              "that step was run. Adds a table slide (not rendered 3D images -- "
                              "see that script's docstring). No-ops if the manifest says skipped.")
    parser.add_argument("--structure-bullets-json", default=None,
                         help="Path to a JSON array of bullet strings about the structure-"
                              "prediction shortlist. Only used if --structure-manifest is given.")
    parser.add_argument("--pptx-name", default="Library_QC_Summary.pptx")
    args = parser.parse_args()

    with open(os.path.join(args.output_dir, "qc_metrics.json")) as f:
        summary_metrics = json.load(f)
    with open(os.path.join(args.output_dir, "plot_manifest.json")) as f:
        manifest = json.load(f)
    with open(args.methodology_json) as f:
        methodology_bullets = json.load(f)

    conclusion_bullets = []
    if args.conclusions_json:
        with open(args.conclusions_json) as f:
            conclusion_bullets = json.load(f)

    ml_images = []
    ml_bullets = []
    if args.ml_manifest:
        with open(args.ml_manifest) as f:
            ml_manifest = json.load(f)
        ml_images = [tuple(x) for x in ml_manifest.get("perplexity_images", [])]
        if args.ml_bullets_json:
            with open(args.ml_bullets_json) as f:
                ml_bullets = json.load(f)

    diversity_estimate = None
    if args.diversity_json:
        with open(args.diversity_json) as f:
            diversity_estimate = json.load(f)

    developability_images = []
    developability_bullets = []
    if args.developability_manifest:
        with open(args.developability_manifest) as f:
            dev_manifest = json.load(f)
        if dev_manifest.get("prevalence_image"):
            developability_images = [["Developability & Liability Screening", dev_manifest["prevalence_image"]]]
        if args.developability_bullets_json:
            with open(args.developability_bullets_json) as f:
                developability_bullets = json.load(f)

    humanness_images = []
    humanness_bullets = []
    if args.humanness_manifest:
        with open(args.humanness_manifest) as f:
            humanness_manifest = json.load(f)
        if not humanness_manifest.get("skipped") and humanness_manifest.get("summary_image"):
            humanness_images = [["Humanness & Immunogenicity Risk (Sapiens)", humanness_manifest["summary_image"]]]
            if args.humanness_bullets_json:
                with open(args.humanness_bullets_json) as f:
                    humanness_bullets = json.load(f)

    structure_table = None
    structure_bullets = []
    if args.structure_manifest:
        with open(args.structure_manifest) as f:
            structure_manifest = json.load(f)
        if not structure_manifest.get("skipped") and structure_manifest.get("rows"):
            # Display file basenames rather than full absolute paths -- a viewer
            # of the deck can't do anything with this machine's local paths
            # anyway; the full paths are still in structure_manifest.json itself.
            headers = structure_manifest["headers"]
            path_col = headers.index("PDB Path") if "PDB Path" in headers else None
            rows = []
            for row in structure_manifest["rows"]:
                row = list(row)
                if path_col is not None:
                    row[path_col] = os.path.basename(row[path_col])
                rows.append(row)
            structure_table = (headers, rows)
            if args.structure_bullets_json:
                with open(args.structure_bullets_json) as f:
                    structure_bullets = json.load(f)

    output_path = os.path.join(args.output_dir, args.pptx_name)
    build_qc_deck(
        output_path=output_path,
        title=args.title,
        subtitle=args.subtitle,
        summary_metrics=summary_metrics,
        methodology_bullets=methodology_bullets,
        distribution_images=[tuple(x) for x in manifest["distribution_images"]],
        pwm_images_by_gene=[tuple(x) for x in manifest["pwm_images"]],
        conclusion_bullets=conclusion_bullets,
        ml_bullets=ml_bullets,
        ml_images=ml_images,
        diversity_estimate=diversity_estimate,
        developability_bullets=developability_bullets,
        developability_images=[tuple(x) for x in developability_images],
        humanness_bullets=humanness_bullets,
        humanness_images=[tuple(x) for x in humanness_images],
        structure_bullets=structure_bullets,
        structure_table=structure_table,
    )
    print(f"Saved {output_path}")


if __name__ == "__main__":
    main()
