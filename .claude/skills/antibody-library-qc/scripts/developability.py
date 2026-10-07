"""
Sequence-based developability / liability screening.

There's no lightweight, pip-installable ML model for this: the real tool
used in the field (TAP, Therapeutic Antibody Profiler) needs a 3D structure
to compute its metrics (CDR charge/hydrophobicity patches on the actual
folded surface), not sequence alone, and isn't pip-installable regardless.
What's here instead is a set of literature-standard sequence motif /
physicochemical heuristic checks -- not predictions from a trained model.
Treat flags as "worth a closer look in the real (structure-aware) tool,"
not verdicts.

Operates on the aligner's already-computed best-alignment window (the
`{Label}_Best_Window` column from qc_pipeline.py / aligner.score_reads_parallel),
so no new alignment work is needed. These are O(design_len) regex/windowed
checks per read -- much cheaper than a model forward pass or the sliding-
window aligner -- so this runs on the full dataset, not a sample.

Also: this pipeline's aligner assumes fixed-length designs (substitutions
only, no indel modeling), so a read's CDR span is always exactly as long as
its design's -- there's no such thing as a "CDR length outlier" to detect
here. That's why length-based liability checks aren't included.
"""

import re

import pandas as pd

KYTE_DOOLITTLE = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5,
    'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8,
    'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
}

DEAMIDATION_RE = re.compile(r'N[GS]')
ISOMERIZATION_RE = re.compile(r'D[GS]')
N_GLYC_RE = re.compile(r'N[^P][ST]')

HYDROPHOBICITY_WINDOW = 5
HYDROPHOBICITY_THRESHOLD = 2.5  # mean Kyte-Doolittle score per residue in window


def screen_developability_liabilities(window_seq, cdr_positions):
    """Flag common sequence-level developability liabilities in one read's
    aligned window. `cdr_positions` is the set of 0-indexed designed/CDR
    positions for this read's assigned design (from a design's
    `allowed_mutations.keys()`, same source plotting.py uses). Returns a
    dict of flags/counts, or None if `window_seq` isn't a usable string."""
    if not isinstance(window_seq, str) or not window_seq:
        return None

    is_cdr = [p in cdr_positions for p in range(len(window_seq))]

    deamidation_total = len(DEAMIDATION_RE.findall(window_seq))
    deamidation_cdr = sum(1 for m in DEAMIDATION_RE.finditer(window_seq) if is_cdr[m.start()])

    isomerization_total = len(ISOMERIZATION_RE.findall(window_seq))
    isomerization_cdr = sum(1 for m in ISOMERIZATION_RE.finditer(window_seq) if is_cdr[m.start()])

    nglyc_total = len(N_GLYC_RE.findall(window_seq))
    nglyc_cdr = sum(1 for m in N_GLYC_RE.finditer(window_seq) if is_cdr[m.start()])

    oxidation_cdr = sum(1 for i, aa in enumerate(window_seq) if is_cdr[i] and aa in ('M', 'W'))
    # A single Met/Trp in a CDR is unremarkable -- natural antibodies have them routinely.
    # Require clustering (2+) before treating it as a flag worth a closer look, not a
    # bare count, or this dwarfs every other check (benchmarked ~80% "flagged" on real
    # data from this alone, which just reflects how common one CDR Met/Trp is).
    oxidation_risk = oxidation_cdr >= 2

    cys_count = window_seq.count('C')
    unpaired_cysteine_risk = (cys_count % 2 != 0)

    hydrophobic_patch = False
    cdr_indices = [i for i, flag in enumerate(is_cdr) if flag]
    for i in range(len(cdr_indices) - HYDROPHOBICITY_WINDOW + 1):
        span = cdr_indices[i:i + HYDROPHOBICITY_WINDOW]
        if span[-1] - span[0] != HYDROPHOBICITY_WINDOW - 1:
            continue  # not contiguous -- span crosses a gap between CDR segments
        mean_kd = sum(KYTE_DOOLITTLE.get(window_seq[p], 0.0) for p in span) / HYDROPHOBICITY_WINDOW
        if mean_kd >= HYDROPHOBICITY_THRESHOLD:
            hydrophobic_patch = True
            break

    deamidation_risk = deamidation_cdr > 0
    isomerization_risk = isomerization_cdr > 0
    nglyc_risk = nglyc_cdr > 0

    any_flag = any([
        deamidation_risk, isomerization_risk, nglyc_risk,
        oxidation_risk, unpaired_cysteine_risk, hydrophobic_patch,
    ])

    return {
        'Deamidation_Motifs_Total': deamidation_total,
        'Deamidation_Motifs_CDR': deamidation_cdr,
        'Deamidation_Risk': deamidation_risk,
        'Isomerization_Motifs_Total': isomerization_total,
        'Isomerization_Motifs_CDR': isomerization_cdr,
        'Isomerization_Risk': isomerization_risk,
        'N_Glyc_Motifs_Total': nglyc_total,
        'N_Glyc_Motifs_CDR': nglyc_cdr,
        'N_Glyc_Risk': nglyc_risk,
        'Oxidation_Prone_CDR_Residues': oxidation_cdr,
        'Oxidation_Risk': oxidation_risk,
        'Cysteine_Count': cys_count,
        'Unpaired_Cysteine_Risk': unpaired_cysteine_risk,
        'CDR_Hydrophobic_Patch': hydrophobic_patch,
        'Any_Liability_Flag': any_flag,
    }


def score_developability_batch(df, design_col, window_col, rules_dict, chain_name):
    """Apply screen_developability_liabilities to every read, using each
    read's already-assigned design's CDR positions. Returns a DataFrame of
    per-read flags (columns prefixed with `chain_name`), aligned to df's
    index -- rows with no usable alignment window get all-NaN flags."""
    cdr_positions_by_gene = {
        gene: set(rules.get('allowed_mutations', {}).keys())
        for gene, rules in rules_dict.items()
    }

    records = []
    for gene, window in zip(df[design_col], df[window_col]):
        cdr_positions = cdr_positions_by_gene.get(gene, set())
        flags = screen_developability_liabilities(window, cdr_positions)
        records.append(flags if flags is not None else {})

    flags_df = pd.DataFrame(records, index=df.index)
    flags_df.columns = [f"{chain_name}_{c}" for c in flags_df.columns]
    return flags_df


def summarize_liability_prevalence(flags_df, read_count_series, chain_name):
    """Abundance-weighted % of reads with each boolean liability flag set,
    among reads that had a usable alignment window. Returns {label: pct}."""
    bool_cols = [c for c in flags_df.columns if c.endswith('_Risk') or c.endswith('_Patch')
                 or c == f"{chain_name}_Any_Liability_Flag"]
    usable = flags_df.notna().all(axis=1)
    total_abundance = read_count_series[usable].sum()
    if total_abundance == 0:
        return {}

    prevalence = {}
    for col in bool_cols:
        flagged_abundance = read_count_series[flags_df[col] == True].sum()
        label = col.replace(f"{chain_name}_", "").replace("_", " ")
        prevalence[label] = round(100 * float(flagged_abundance) / float(total_abundance), 2)
    return prevalence


def plot_liability_prevalence(prevalence_by_chain, output_path):
    """Grouped bar chart: liability flag (x) vs. abundance-weighted
    prevalence % (y), one bar group per chain in `prevalence_by_chain`
    (e.g. {"Heavy": {...}, "Light": {...}})."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    labels = sorted({label for prevalence in prevalence_by_chain.values() for label in prevalence})
    chains = list(prevalence_by_chain.keys())
    x = np.arange(len(labels))
    width = 0.8 / max(len(chains), 1)

    fig, ax = plt.subplots(figsize=(9, 4.5))
    for i, chain in enumerate(chains):
        values = [prevalence_by_chain[chain].get(label, 0.0) for label in labels]
        ax.bar(x + i * width - (len(chains) - 1) * width / 2, values, width, label=chain)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha='right', fontsize=8)
    ax.set_ylabel('Abundance-Weighted Prevalence (%)', fontsize=10)
    ax.set_title('Sequence-Based Developability Liability Flags', fontsize=11, fontweight='bold')
    ax.legend(title='Chain')
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return output_path
