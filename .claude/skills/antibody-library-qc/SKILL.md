---
name: antibody-library-qc
description: QC and analyze antibody (or nanobody/scFv) display library NGS sequencing data — scoring reads against a design-rule workbook for intended CDR/designed mutations vs true framework errors, computing valid-read and paired-pairing efficiency, flagging synthesis/cloning error hotspots, quantifying design bias via KL divergence, scoring sequence fitness with an antibody-specific protein language model (AbLang2), screening sequence-based developability/liability flags, scoring humanness/immunogenicity risk (Sapiens), predicting 3D structure for a shortlist (IgFold), estimating true library diversity (Chao1 / Good-Turing) from raw read counts, and assembling a PowerPoint QC report. Use this whenever the user has a CSV/TSV of sequenced reads (amino acid sequences + read counts), optionally with a design workbook of wild-type sequences and allowed substitutions, and wants library QC metrics, error-profile histograms, framework-vs-CDR PWM plots, design-bias/fitness scoring, developability/liability screening, humanness/immunogenicity scoring, structure prediction, diversity/complexity estimates, or a QC summary deck — even if they just say "check my library sequencing data," "how clean is this library," "score these reads against the design," "how much diversity does this library actually have," "is my library biased," "will this immunogenic," "any developability red flags," or hand you NGS output without using any of these words explicitly.
---

# Antibody Library QC

Generalized from three working notebook pipelines in this repo
(`notebooks/02_library_analysis_v2.ipynb`, `notebooks/03_ml_enrichment_and_bias.ipynb`,
`notebooks/NGS_Library_Diversity_Estimation.ipynb`) that QC'd a paired heavy/light
antibody library, then extended with three more ML-driven analyses (developability
screening, humanness scoring, structure prediction) added after the fact. The
approach applies to any display-library NGS QC job: single-domain (nanobody/VHH),
scFv, or paired Fab/IgG libraries, as long as you have (1) observed reads and (2) a
wild-type + allowed-mutation design spec to score them against. There are six
independent analyses below — run whichever ones answer the question actually being
asked, not all six by default.

## When to use this

Reach for this when the user wants to know **how clean a sequenced library is relative
to what was designed**, **whether its composition drifted from what was intended**, or
**how much true diversity it actually contains** — not just arbitrary sequence analysis.
Signals: they mention read counts/abundance, a design or PWM spreadsheet, CDR vs
framework errors, synthesis or cloning error rates, "did the library build correctly,"
design bias, sequence fitness/stability/naturalness, or how many unique
clones/variants a library really has versus how many were actually sequenced.

Six independent analyses, pick based on the question:

| Question | Analysis | Needs a design workbook? | Extra setup? |
|---|---|---|---|
| "How clean is this library / did it build correctly?" | 1. Core QC | Yes | No |
| "Is the library's amino acid composition biased vs. what we designed? Do framework errors hurt predicted fitness?" | 2. Design bias & ML fitness | Yes, plus Core QC's scored output | No |
| "How many unique clones does this library really have?" | 3. Diversity estimation | No — raw read counts only | No |
| "Any developability red flags (deamidation, oxidation, hydrophobic patches, ...)?" | 4. Developability screening | Yes, plus Core QC's scored output | No |
| "Will this be immunogenic / how human does it look?" | 5. Humanness scoring (Sapiens) | Yes, plus Core QC's scored output | Yes — separate Python 3.8 env |
| "What does this antibody's 3D structure look like?" | 6. Structure prediction (IgFold) | Yes, plus Core QC's scored output | Yes — separate Python 3.11 env |

## The core idea

Observed reads rarely align exactly to the designed window — there's untrimmed primer,
leader sequence, or constant-domain carryover at the ends. The pipeline's aligner
(`scripts/aligner.py`) slides the full-length wild-type design across each read and
keeps the best-fitting window, then splits mismatches in that window into two buckets:

- **Intended mutations** — mismatches at a designed position, in that position's
  allowed-amino-acid set (i.e. the library did what it was designed to do).
- **True errors** — everything else (synthesis errors, cloning errors, PCR/primer
  scarring, sequencing noise).

Everything downstream (efficiency metrics, error histograms, hotspot flags) is built on
this read-level classification.

## Design workbook format

`scripts/pwm_extraction.py` parses a specific but common spreadsheet convention (matches
this repo's `AffMat_design.xlsx`): one block per gene, stacked vertically in a sheet,
where a `WT` marker cell flags the row/column where each designed position's allowed
amino acids and frequencies are tabulated. Read the docstring at the top of that file
before assuming a new workbook matches — if the layout differs, adjust the row/column
offsets there rather than rewriting the call sites in `qc_pipeline.py`.

If the user's design info isn't in this spreadsheet shape at all (e.g. a plain list of
allowed mutations per position), skip `pwm_extraction.py` and construct the
`{design_name: {'wt_sequence': ..., 'allowed_mutations': {pos: [aa, ...]}}}` dict
directly — that's the only interface `aligner.py` and `plotting.py` actually need.

## Analysis 1: Core QC

1. **Run the scoring pipeline.** This does all the deterministic data crunching —
   nothing here requires interpretation, so just run it:

   ```bash
   python scripts/qc_pipeline.py \
     --sequences-csv <reads.csv> \
     --design-xlsx <design.xlsx> \
     --heavy-sheet "Heavy Chain Design" --heavy-col HAA \
     --light-sheet "Light Chain Design" --light-col LAA \
     --read-count-col read_count \
     --error-threshold 2 \
     --fidelity-threshold 0.85 \
     --output-dir <output_dir>
   ```

   Single-chain libraries (nanobody, scFv) just omit the `--light-*` flags and pass
   whatever the one chain's sheet/column names are via `--heavy-sheet`/`--heavy-col`
   (the "heavy" label is just the first-chain slot, not a biology assumption).

   `--error-threshold` is the max true framework errors a read can carry and still
   count as "valid" — 2 is a reasonable default for most NGS error rates but ask the
   user if they have a different QC bar. `--fidelity-threshold` controls when a
   framework position gets flagged as a hotspot in the PWM plots (default: WT
   observed in <85% of reads at that position).

   This writes `qc_metrics.json`, `plot_manifest.json`, `scored_reads.csv`, the
   abundance-weighted error-distribution histograms, and per-gene split
   Framework-vs-CDR PWM plots into `--output-dir`.

2. **Actually look at the results before writing anything about them.** Read
   `qc_metrics.json` and view the generated PNGs (histograms + split PWM plots).
   Efficiency numbers below ~50-60% valid reads, large gaps between heavy/light
   individual efficiency and paired efficiency, or a visible secondary peak in the
   true-error histogram are all worth calling out. Framework positions with red
   asterisks in the split PWM plots are literal error hotspots — look at which FW
   segment they fall in and whether they cluster at one terminus (classic sign of
   primer mis-synthesis or scarring) versus scattered randomly (sign of general
   synthesis noise).

3. **Report findings to the user directly**, or if they want a shareable deck, write
   two small JSON files — a methodology bullet list and a conclusions bullet list
   containing only things actually observed in *this* run's metrics/plots, not
   boilerplate carried over from a different dataset — then build the deck:

   ```bash
   python scripts/build_report.py \
     --output-dir <output_dir> \
     --title "Antibody Library Analysis" \
     --subtitle "Framework & CDR Split PWM Profiling" \
     --methodology-json methodology.json \
     --conclusions-json conclusions.json \
     --pptx-name Library_QC_Summary.pptx
   ```

   Never invent conclusion bullets before step 2 has actually happened — the whole
   point of this report is that the numbers in it are real.

## Analysis 2: Design bias & ML fitness scoring

Builds on Analysis 1's output — run `qc_pipeline.py` first so `scored_reads.csv` exists.
This answers two different questions than the core QC does: not "is this read within
error tolerance" but "did the library's *overall composition* drift from what was
designed" and "do framework errors actually correlate with worse predicted sequence
fitness, or are they just noise."

```bash
python scripts/ml_bias_fitness.py \
  --scored-reads-csv <output_dir>/scored_reads.csv \
  --design-xlsx <design.xlsx> \
  --heavy-sheet "Heavy Chain Design" --heavy-col HAA \
  --light-sheet "Light Chain Design" --light-col LAA \
  --read-count-col read_count \
  --output-dir <output_dir>
```

Both of the following run **by default** — no extra flags needed:

- **KL divergence** (cheap, no extra dependencies beyond `scipy`): for every
  designed/CDR position, compares the amino acid frequencies the design *intended*
  (from the PWM) against what sequencing *actually observed* at that position. High
  divergence flags positions where the real library composition diverged from the
  design — commonly oligo/trimer synthesis bias toward certain codons. Written to
  `kl_divergence.csv`; the top rows by `KL_Divergence` are the ones worth mentioning.
  Always runs on the full dataset, not a sample.
- **Fitness-LM pseudo-perplexity**: scores each sequence with a protein language
  model. Lower perplexity = the model finds the sequence more natural/stable. The
  scatter plot this produces (true framework errors vs. perplexity) tells you
  whether framework errors are actually destabilizing or just harmless noise.
  Scores the **aligned `{Label}_Best_Window`**, not the raw sequence column — raw
  reads carry CDR3/J-segment/constant-domain carryover that the design's
  `wt_sequence` (used for the fold-compatible reference line) never includes, so
  scoring the untrimmed read against a baseline built from the trimmed WT would
  inflate every read's perplexity relative to WT regardless of actual error count.
  (This was a real bug here briefly — caught by a user noticing that even
  zero-error reads scored far above the WT line, which shouldn't happen once the
  comparison is apples-to-apples. `humanness_scoring.py`/Analysis 5 has the same
  fix for the same reason.)
  Defaults to **AbLang2** (`--fitness-backend ablang2`), an antibody-specific LM
  (trained on real antibody repertoires, not all of UniRef the way ESM-2 is) —
  pass `--fitness-backend esm2` to use the earlier general-purpose model instead,
  kept around for comparison. AbLang2 scores heavy+light **jointly** by design; this
  script uses it to score each chain alone too (so the existing per-chain plots keep
  working) and *also* adds a third plot — `Paired_Fitness_Perplexity`, a joint
  "pairing compatibility" score with its own WT-pair baseline — only produced when
  both chains are given, since that's specifically about how well the two chains fit
  *together*, not just individually.

  Because this is a transformer forward pass per sequence — orders of magnitude
  slower than the alignment work — **`ml_bias_fitness.py` automatically caps it to a
  random sample whenever the dataset is bigger than that** (`--fitness-sample-size`;
  defaults to 15,000 for ablang2 or 50,000 for esm2, reflecting AbLang2's ~5-6x
  higher per-call cost and the extra joint-pairing call). The sample is seeded
  (`--fitness-sample-seed`, default 42) so repeat runs are reproducible. Pass
  `--skip-fitness-lm` to disable this step entirely, or `--fitness-sample-size 0` to
  force scoring every row regardless of size (expect hours past roughly 100k rows on
  CPU, for either backend). `ml_manifest.json` records whether sampling kicked in
  (`fitness_backend`, `fitness_sampled`, `fitness_rows_scored`, `fitness_total_rows`)
  so you can report that honestly rather than implying every row was scored.

Always pass the resulting manifest into `build_report.py` (next section) so the KL
divergence findings and the fitness-LM scatter plots land in the deck automatically —
don't treat Analysis 2 as a one-off extra that has to be requested separately; running
it and folding it into the deck is the default end-to-end flow now. Pass the manifest
it writes (plus your own bullets summarizing what the KL/perplexity results actually
showed) to `build_report.py`:

```bash
python scripts/build_report.py \
  --output-dir <output_dir> \
  --methodology-json methodology.json \
  --ml-manifest <output_dir>/ml_manifest.json \
  --ml-bullets-json ml_bullets.json \
  --pptx-name Library_QC_Summary.pptx
```

## Analysis 3: Library diversity estimation

Standalone — no design workbook, no dependency on Analysis 1 or 2. Answers "how many
distinct clones does this library actually contain, including ones sequencing missed
entirely?" using two nonparametric richness estimators (Chao1, and a Good-Turing-style
Poisson coverage model) built from just the frequency-of-frequencies of read counts
(mainly singletons `f1` and doubletons `f2`). Use this on the raw/unfiltered sequencing
output (every observed unique clone with its read count), not a design-rule-filtered
subset — the estimators assume the input reflects the full observed clone population.

```bash
python scripts/diversity_estimation.py \
  --input-csv <reads.csv or reads.csv.gz> \
  --read-count-col read_count \
  --output-dir <output_dir>
```

Streams the file in chunks (`--chunksize`, default 100,000 rows) so multi-million-row
gzipped files work without loading everything into memory. Writes
`diversity_estimate.json` with both point estimates and confidence intervals at 80/90/95/99%.

A library where singletons vastly outnumber doubletons (`f1 >> f2`) means sequencing
barely sampled the true diversity — that's exactly when Chao1/Poisson estimates will
diverge furthest from the naive "unique rows observed" count, and worth flagging to
the user rather than reporting the observed count as if it were the true diversity.

To add this to the deck from Analysis 1, pass its output straight through to
`build_report.py` — it doesn't need to live under the same `--output-dir` since this
analysis is independent of the design-rule scoring. Rather than its own slide buried
later in the deck, the point estimates and confidence-interval table are folded onto
the front "Library Build Summary & Pairing Efficiency" slide, so every top-line metric
— QC efficiency and diversity alike — is visible together at the start of the deck:

```bash
python scripts/build_report.py \
  --output-dir <output_dir> \
  --methodology-json methodology.json \
  --diversity-json <diversity_output_dir>/diversity_estimate.json \
  --pptx-name Library_QC_Summary.pptx
```

## Analysis 4: Developability & liability screening

Builds on Analysis 1's output. Answers "does this sequence have known red-flag motifs
worth a closer look" — deamidation (`N[GS]`), isomerization (`D[GS]`), oxidation-prone
CDR residue clustering (2+ Met/Trp in a CDR — a single one is unremarkable and *not*
flagged, see `developability.py`'s comment on why that threshold matters), an
N-glycosylation motif (`N[^P][ST]`), cysteine-count anomalies, and CDR hydrophobic
patches (Kyte-Doolittle windowed score).

**This is not a trained ML model** — there's no lightweight pip-installable one, since
the real tool (TAP, Therapeutic Antibody Profiler) needs an actual 3D structure to
compute its metrics (real surface patches), not sequence alone. These are literature-
standard sequence/physicochemical heuristics instead. Treat flags as "worth checking
in a real structure-aware tool," not verdicts — see `developability.py`'s module
docstring for the full caveat.

```bash
python scripts/developability_screen.py \
  --scored-reads-csv <output_dir>/scored_reads.csv \
  --design-xlsx <design.xlsx> \
  --heavy-sheet "Heavy Chain Design" --heavy-col HAA \
  --light-sheet "Light Chain Design" --light-col LAA \
  --read-count-col read_count \
  --output-dir <output_dir>
```

No model, no re-alignment (reuses the cached `{Label}_Best_Window` column) — cheap
enough to run on the **full** dataset, not a sample. Writes `developability_flags.csv`
(per-read), `developability_summary.json` (abundance-weighted prevalence % per flag),
and a prevalence bar chart. Fold into the deck the same way as the others:

```bash
python scripts/build_report.py \
  --output-dir <output_dir> \
  --methodology-json methodology.json \
  --developability-manifest <output_dir>/developability_manifest.json \
  --developability-bullets-json developability_bullets.json \
  --pptx-name Library_QC_Summary.pptx
```

## Environment setup for Analyses 5 & 6

Sapiens (humanness) and IgFold (structure) both need Python versions incompatible
with this skill's main venv (3.9.6) *and* with each other — Sapiens needs <=3.8,
IgFold needs >=3.11 — so each lives in its own dedicated environment that the main
venv calls into via `subprocess`. **Both bridge scripts skip gracefully with a clear
message if their environment isn't set up** — this is a one-time, manual, per-machine
prerequisite, not something every run can assume. See
`references/environment_setup.md` for the exact bootstrap commands, including the
workaround needed on this machine for Sapiens' Python 3.8 (building from source
segfaults on modern macOS/Apple Silicon; a prebuilt CPython distribution works). Both
venv paths are auto-detected relative to the repo root
(`venv-sapiens38/bin/python3.8`, `venv-igfold311/bin/python`) — override with
`--sapiens-python`/`--igfold-python` if yours live elsewhere. Both directories are
already in `.gitignore`.

## Analysis 5: Humanness / immunogenicity risk (Sapiens)

Builds on Analysis 1's output. Answers "will this look foreign to the immune system" —
scores sequences with Sapiens, a BERT-style LM trained specifically on natural human
antibody repertoires, deriving a perplexity-style score (lower = more human-like,
same framing as the fitness-LM metrics) and a simpler "% positions matching the
model's top-1 predicted human residue," closer to how OASis/Sapiens humanness is
usually reported in the literature. Also scores each design's own WT sequence as a
baseline, same idea as the fold-compatibility reference line elsewhere in this skill.

```bash
python scripts/humanness_scoring.py \
  --scored-reads-csv <output_dir>/scored_reads.csv \
  --design-xlsx <design.xlsx> \
  --heavy-sheet "Heavy Chain Design" --heavy-col HAA \
  --light-sheet "Light Chain Design" --light-col LAA \
  --output-dir <output_dir>
```

Same cost reasoning and auto-sampling as the fitness LM (`--sample-size`, default
15,000; `--sample-seed`, default 42). Writes `humanness_summary.json` and a summary
bar chart (sample mean vs. WT baseline, both metrics, both chains). Fold into the deck:

```bash
python scripts/build_report.py \
  --output-dir <output_dir> \
  --methodology-json methodology.json \
  --humanness-manifest <output_dir>/humanness_manifest.json \
  --humanness-bullets-json humanness_bullets.json \
  --pptx-name Library_QC_Summary.pptx
```

## Analysis 6: Antibody structure prediction (IgFold)

Builds on Analysis 1's output. Predicts an actual 3D structure (not just a score) for
a **small, explicit shortlist** — default the top 25 most-abundant pairs that pass
both chains' error threshold — using IgFold, an antibody-specific structure predictor.
Never run this on the full dataset or a large sample: even though it benchmarked fast
here (~1.2s/structure on CPU, no refinement), structure prediction is categorically
heavier than everything else in this skill, and nothing about running it at NGS scale
has been validated.

```bash
python scripts/structure_prediction.py \
  --scored-reads-csv <output_dir>/scored_reads.csv \
  --heavy-col HAA --light-col LAA --read-count-col read_count \
  --heavy-error-col Heavy_True_Errors --light-error-col Light_True_Errors \
  --shortlist-size 25 \
  --output-dir <output_dir>
```

There's no lightweight way to render a 3D structure into a static PowerPoint image, so
results land in the deck as a **table** (sequence ID, abundance, mean predicted RMSD
in Å — IgFold's own per-structure confidence estimate, lower = more confident — and
the PDB file basename), not a picture. The actual PDB files are written to
`<output_dir>/structures/` for anyone who wants to open them in a real structure
viewer. IgFold's code/weights are licensed for **non-commercial use only** (JHU
Academic Software License) — flag this to the user if commercial use is in play; see
`igfold_worker.py`'s docstring for the license link.

```bash
python scripts/build_report.py \
  --output-dir <output_dir> \
  --methodology-json methodology.json \
  --structure-manifest <output_dir>/structure_manifest.json \
  --structure-bullets-json structure_bullets.json \
  --pptx-name Library_QC_Summary.pptx
```

## Script reference

- `scripts/pwm_extraction.py` — `extract_full_pwm(filepath, sheet_name)`, parses the
  design workbook into `{gene: {wt_sequence, pwm, allowed_mutations}}`.
- `scripts/aligner.py` — `score_read_against_designs(read_seq, rules_dict)`, the
  sliding-window classifier described above, plus `score_reads_parallel(...)` which
  fans scoring out across CPU cores. This is the single most expensive step in the
  whole skill (it tests every read against every design), so for large files always
  go through `score_reads_parallel` rather than calling `score_read_against_designs`
  in a plain `.apply()` loop. See "Performance at scale" below.
- `scripts/plotting.py` — `plot_error_distributions(...)` for the abundance-weighted
  histograms, `generate_split_observed_pwm(...)` for the per-gene FW-vs-CDR charts.
- `scripts/qc_pipeline.py` — CLI that wires the above together, computes efficiency
  metrics, and writes `qc_metrics.json` / `plot_manifest.json` / `scored_reads.csv`.
- `scripts/pptx_report.py` — `build_qc_deck(...)`, pure presentation assembly from
  already-computed metrics/images/bullets (including the optional ML slide/images).
- `scripts/build_report.py` — CLI wrapper around `build_qc_deck` that reads
  `qc_metrics.json`/`plot_manifest.json`, the methodology/conclusions JSON files, and
  optionally `ml_manifest.json` from `ml_bias_fitness.py`.
- `scripts/bias_and_fitness.py` — `calculate_pwm_kl_divergence(...)` plus both
  fitness-LM backends: `load_esm2`/`compute_esm2_perplexity`/`score_esm2_perplexity_batch`
  and `load_ablang2`/`compute_ablang2_perplexity`/`score_ablang2_chain_batch`/
  `score_ablang2_paired_batch`, each with a matching `compute_*_wt_reference(_perplexity)`
  baseline helper, plus the shared `plot_perplexity_vs_errors(...)`. The Analysis 2 logic,
  reusable outside the CLI.
- `scripts/ml_bias_fitness.py` — CLI that runs Analysis 2 end to end (KL divergence +
  auto-sampled fitness-LM scoring, both on by default — see Analysis 2 above) and
  writes `kl_divergence.csv` / `ml_manifest.json` / `scored_reads_with_fitness.csv`
  (or `..._sample.csv` if sampling kicked in).
- `scripts/diversity_estimation.py` — `estimate_diversity(...)`, `chao1_multi_ci(...)`,
  `poisson_solver(...)`, `poisson_multi_ci(...)`, plus a CLI (Analysis 3). Fully
  standalone, no imports from the rest of this skill.
- `scripts/developability.py` — `screen_developability_liabilities(...)`,
  `score_developability_batch(...)`, `summarize_liability_prevalence(...)`,
  `plot_liability_prevalence(...)`: the Analysis 4 logic, reusable outside the CLI.
- `scripts/developability_screen.py` — CLI that runs Analysis 4 end to end and writes
  `developability_flags.csv` / `developability_summary.json` /
  `developability_manifest.json`.
- `scripts/sapiens_worker.py` — standalone, runs **in the Sapiens Python 3.8
  environment only** (not importable from the main venv). Scores sequences with
  Sapiens; invoked via subprocess by `humanness_scoring.py`.
- `scripts/humanness_scoring.py` — CLI that runs **in the main venv**, orchestrates
  Analysis 5 end to end (auto-sampling, subprocess bridge to the Sapiens env, WT
  baseline, summary chart) and writes `humanness_summary.json` /
  `humanness_manifest.json`. Skips gracefully if the Sapiens env isn't set up.
- `scripts/igfold_worker.py` — standalone, runs **in the IgFold Python 3.11
  environment only**. Folds one structure per shortlisted pair; invoked via
  subprocess by `structure_prediction.py`.
- `scripts/structure_prediction.py` — CLI that runs **in the main venv**, orchestrates
  Analysis 6 end to end (shortlist selection, subprocess bridge to the IgFold env) and
  writes `structure_manifest.json`. Skips gracefully if the IgFold env isn't set up.

All scripts were validated against this repo's `preview_1000_pairs.csv` +
`AffMat_design.xlsx`. Analysis 1 (291 valid heavy reads, 326 valid light, 102
fully-paired, 101 unique quality-filtered pairs) and Analysis 3 (same Chao1/Poisson
estimates given the same input) reproduce their source notebooks' numbers exactly.
Analysis 2 (KL divergence, fitness LM) is *intentionally* slightly different from
`03_ml_enrichment_and_bias.ipynb` when run as part of this pipeline — see
"Performance at scale" below for why, and why that's a correctness improvement, not
a regression. Analyses 4-6 are new (not from a source notebook) and were smoke-tested
end to end on the preview dataset, including both subprocess bridges.

## Performance at scale

`aligner.py`'s sliding-window search is the dominant cost in every analysis here: it
tests every read against every design. On a benchmark against this repo's real
~3.7M-row `clean_pairs_annotated.csv.gz`, single-process scoring ran at roughly
1-1.5ms/read for one chain — call it 1-2 hours per chain on a single core, with
`plotting.py` and `bias_and_fitness.py` each separately re-deriving the same
alignment on top of that if run the old way. Two things make this tractable on
large files:

1. **`score_reads_parallel` instead of `.apply()`.** Each read is scored
   independently, so this is embarrassingly parallel. `qc_pipeline.py`'s
   `score_chain()` always goes through `score_reads_parallel`, which spreads reads
   across `--n-processes` cores (default: all of them) via `multiprocessing.Pool`,
   giving close to linear speedup with core count. For a few hundred reads or fewer
   it just runs in-process — not worth paying process-startup cost at that size.

2. **Cache the winning alignment window, reuse it everywhere.** `score_reads_parallel`
   doesn't just return which design/how many errors — it returns the actual
   best-matching substring of the read (`{Label}_Best_Window`, written into
   `scored_reads.csv`). `plotting.py`'s `generate_split_observed_pwm` and
   `bias_and_fitness.py`'s `calculate_pwm_kl_divergence` both accept a
   `best_window_col` argument: when it's given (which `qc_pipeline.py` and
   `ml_bias_fitness.py` do automatically) they reuse that cached window instead of
   re-running their own sliding-window search from scratch, cutting out a
   meaningful chunk of redundant work on top of the main scoring pass. They still
   fall back to recomputing from the raw sequence column if used standalone without
   that column.

   **Why this changes Analysis 2's numbers slightly**: the aligner picks a window by
   minimizing *true errors only* (an allowed CDR substitution doesn't count against
   a window). The old standalone window search in `plotting.py`/`bias_and_fitness.py`
   (and the original notebooks) instead minimized *total mismatches*, counting
   intended CDR mutations as if they were errors for the purpose of picking the best
   window. Those two criteria agree almost everywhere, but can disagree right at
   designed/CDR positions on reads that actually carry a real mutation there — which
   is exactly where the two sets of outputs differ. Reusing the cached window makes
   every analysis in this skill agree on one definition of "best alignment" (the same
   one the QC error/intended counts are based on), which is more internally
   consistent than the original notebooks' behavior, not less correct — but don't
   expect bit-for-bit identical PWM/KL numbers to `03_ml_enrichment_and_bias.ipynb`
   on data with real CDR diversity.

For files in the preview_1000_pairs.csv size range this is all moot (seconds either
way). It starts to matter past tens of thousands of rows, and matters a lot past
hundreds of thousands.

One more thing worth knowing for real NGS exports: `qc_pipeline.py` only reads the
columns it actually uses (`--heavy-col`/`--light-col`/`--read-count-col`) via
`usecols`, not the whole file. Real annotated exports often carry dozens of
provenance/allele-call columns alongside `HAA`/`LAA`/`read_count` — loading and then
re-writing all of them into `scored_reads.csv` would multiply memory use and disk
size for columns nothing downstream touches.

## Dependencies

**Main venv** (`venv/`, Python 3.9.6 — activate before running any script except the
two standalone workers): Core QC (Analysis 1) and report assembly need `pandas`,
`numpy`, `matplotlib`, `seaborn`, `openpyxl`, `python-pptx`. Analysis 2 adds `scipy`
always, plus `torch`/`ablang2` (default backend) or `torch`/`transformers` (`--fitness-backend
esm2`) — both download model weights on first use, so that step needs network access
the first time it runs. Analysis 4 (developability): no extra dependencies, just
`pandas`/regex. Analyses 5 & 6's *bridge* scripts (`humanness_scoring.py`,
`structure_prediction.py`) also run in the main venv and need nothing beyond
`pandas` — the heavy dependencies live in their own environments. Analysis 3
(diversity): `pandas` only.

**Sapiens environment** (`venv-sapiens38/`, Python 3.8, see "Environment setup"
above): `sapiens`, `torch`, `pandas`. Only `sapiens_worker.py` runs here.

**IgFold environment** (`venv-igfold311/`, Python 3.11, see "Environment setup"
above): `igfold` (pulls in `antiberty`, `torch`, `biopython`), `pandas`. Only
`igfold_worker.py` runs here.
