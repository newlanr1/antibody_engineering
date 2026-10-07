"""
Plotting helpers for antibody library QC: abundance-weighted error-profile
histograms, and per-gene split Framework-vs-CDR observed-PWM charts.
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def plot_error_distributions(df, intended_col, error_col, read_count_col,
                              chain_name, output_path,
                              intended_color='mediumseagreen', error_color='crimson'):
    """Save a 2-panel figure: intended-mutation histogram and true-error
    histogram, both weighted by read_count_col so abundant reads count more."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 4.5))

    intended_max = df[intended_col].max(skipna=True)
    intended_max = 0 if pd.isna(intended_max) else int(intended_max)
    sns.histplot(data=df, x=intended_col, weights=read_count_col,
                 bins=range(0, intended_max + 2), color=intended_color, ax=axes[0])
    axes[0].set_title(f'{chain_name}: Intended Designed Mutations', fontsize=12)
    axes[0].set_xlabel('Number of Successful Mutations per Read')

    error_max = df[error_col].max(skipna=True)
    error_max = 0 if pd.isna(error_max) else int(error_max)
    sns.histplot(data=df, x=error_col, weights=read_count_col,
                 bins=range(0, error_max + 2), color=error_color, ax=axes[1])
    axes[1].set_title(f'{chain_name}: True Error Profile', fontsize=12)
    axes[1].set_xlabel('Number of True Errors per Read')
    axes[1].set_xlim(0, max(20, error_max))

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return output_path


def generate_split_observed_pwm(df, rules_dict, chain_col, design_col, read_count_col,
                                 chain_name, output_dir, fidelity_threshold=0.85,
                                 best_window_col=None):
    """For every gene in rules_dict, compute observed per-position amino acid
    frequencies across all reads assigned to that gene and render a stacked
    2-panel bar chart: Framework regions on top, CDR (designed) regions on
    the bottom. Framework positions below `fidelity_threshold` wild-type
    fidelity are flagged with a red asterisk.

    If `best_window_col` is given (the window column aligner.score_reads_parallel
    already computed, e.g. "Heavy_Best_Window"), reuses that cached alignment
    instead of re-running the sliding-window search here — every read has
    already been aligned to its winning design once, so there's no reason to
    redo it per-gene. Omit it to fall back to recomputing from `chain_col`
    (e.g. if this is being called standalone, without qc_pipeline.py).

    Returns {gene_name: image_path} for genes with at least one scorable read.
    """
    os.makedirs(output_dir, exist_ok=True)
    image_paths = {}

    for gene_name, rules in rules_dict.items():
        wt_seq = rules['wt_sequence']
        design_len = len(wt_seq)
        allowed_muts = rules.get('allowed_mutations', {})

        cdr_positions = set(p for p in allowed_muts.keys() if p < design_len)
        fw_positions = set(range(design_len)) - cdr_positions

        gene_reads = df[df[design_col] == gene_name]
        if len(gene_reads) == 0:
            continue

        aa_list = sorted(list("ACDEFGHIKLMNPQRSTVWY"))
        pos_counts = {p: {aa: 0 for aa in aa_list} for p in range(design_len)}
        total_depth = 0

        for _, row in gene_reads.iterrows():
            read_count = row[read_count_col]
            best_window = row[best_window_col] if best_window_col else None

            if not isinstance(best_window, str):
                read_seq = row[chain_col]
                if not isinstance(read_seq, str) or len(read_seq) < design_len:
                    continue

                best_window = None
                min_err = float('inf')
                for i in range(len(read_seq) - design_len + 1):
                    window = read_seq[i:i + design_len]
                    errs = sum(1 for pos in range(design_len) if window[pos] != wt_seq[pos])
                    if errs < min_err:
                        min_err = errs
                        best_window = window

            if best_window:
                total_depth += read_count
                for pos in range(design_len):
                    res = best_window[pos]
                    if res in pos_counts[pos]:
                        pos_counts[pos][res] += read_count

        if total_depth == 0:
            continue

        freq_df = pd.DataFrame(pos_counts).T / total_depth
        sorted_fw = sorted(list(fw_positions))
        sorted_cdr = sorted(list(cdr_positions))
        colors = plt.cm.tab20(np.linspace(0, 1, 20))

        fig, (ax_fw, ax_cdr) = plt.subplots(2, 1, figsize=(11, 7.2))

        # --- TOP PLOT: Framework Regions (FW) ---
        bottom_fw = np.zeros(len(sorted_fw))

        for idx, aa in enumerate(aa_list):
            values = [freq_df.loc[p, aa] if p in freq_df.index else 0.0 for p in sorted_fw]
            ax_fw.bar(range(len(sorted_fw)), values, bottom=bottom_fw, label=aa,
                      color=colors[idx % 20], width=0.75)
            bottom_fw += np.array(values)

        for i, pos in enumerate(sorted_fw):
            if pos in freq_df.index:
                wt_aa = wt_seq[pos]
                wt_freq = freq_df.loc[pos, wt_aa] if wt_aa in freq_df.columns else 0
                if wt_freq < fidelity_threshold:
                    ax_fw.text(i, 1.02, '*', ha='center', va='bottom', color='crimson',
                               fontsize=11, fontweight='bold')

        fw_segment_idx = 1
        segment_start = 0
        for i in range(1, len(sorted_fw)):
            if sorted_fw[i] != sorted_fw[i - 1] + 1:
                ax_fw.axvline(i - 0.5, color='black', linestyle='--', linewidth=0.8, alpha=0.7)
                mid_point = (segment_start + i - 1) / 2
                ax_fw.text(mid_point, 1.08, f"FW{fw_segment_idx}", ha='center', va='bottom',
                           fontsize=8, fontweight='bold', color='navy')
                fw_segment_idx += 1
                segment_start = i
        if sorted_fw:
            mid_point = (segment_start + len(sorted_fw) - 1) / 2
            ax_fw.text(mid_point, 1.08, f"FW{fw_segment_idx}", ha='center', va='bottom',
                       fontsize=8, fontweight='bold', color='navy')

        ax_fw.set_title(f"{chain_name} Design: {gene_name} — Framework Regions", fontsize=11,
                        fontweight='bold', pad=25)
        ax_fw.set_ylabel('Observed Freq', fontsize=9)
        ax_fw.set_ylim(0, 1.22)
        ax_fw.set_xticks(range(len(sorted_fw)))
        ax_fw.set_xticklabels([wt_seq[p] if p < len(wt_seq) else '-' for p in sorted_fw], fontsize=6.5)
        ax_fw.set_xlabel('Wild-Type Amino Acid Residue', fontsize=8)

        ax_fw_top = ax_fw.twiny()
        ax_fw_top.set_xlim(ax_fw.get_xlim())
        ax_fw_top.set_xticks(range(len(sorted_fw)))
        ax_fw_top.set_xticklabels([f"{p + 1}" for p in sorted_fw], fontsize=5.5, rotation=90)
        ax_fw_top.set_xlabel('Sequence Position (1-indexed)', fontsize=8, labelpad=5)

        # --- BOTTOM PLOT: Designed (e.g. CDR) Regions ---
        bottom_cdr = np.zeros(len(sorted_cdr))
        cdr_labels = [f"P{p + 1}\n({wt_seq[p]})" if p < len(wt_seq) else f"P{p + 1}\n(-)" for p in sorted_cdr]

        for idx, aa in enumerate(aa_list):
            values = [freq_df.loc[p, aa] if p in freq_df.index else 0.0 for p in sorted_cdr]
            ax_cdr.bar(range(len(sorted_cdr)), values, bottom=bottom_cdr, label=aa,
                       color=colors[idx % 20], width=0.75)
            bottom_cdr += np.array(values)

        cdr_segment_idx = 1
        cdr_start = 0
        for i in range(1, len(sorted_cdr)):
            if sorted_cdr[i] != sorted_cdr[i - 1] + 1:
                ax_cdr.axvline(i - 0.5, color='black', linestyle='--', linewidth=0.8, alpha=0.7)
                mid_point = (cdr_start + i - 1) / 2
                ax_cdr.text(mid_point, 1.03, f"CDR{cdr_segment_idx}", ha='center', va='bottom',
                            fontsize=8, fontweight='bold', color='darkgreen')
                cdr_segment_idx += 1
                cdr_start = i
        if sorted_cdr:
            mid_point = (cdr_start + len(sorted_cdr) - 1) / 2
            ax_cdr.text(mid_point, 1.03, f"CDR{cdr_segment_idx}", ha='center', va='bottom',
                        fontsize=8, fontweight='bold', color='darkgreen')

        ax_cdr.set_title(f"{chain_name} Design: {gene_name} — Designed (CDR) Regions", fontsize=11,
                         fontweight='bold')
        ax_cdr.set_xlabel('Position & Wild-Type Residue', fontsize=9)
        ax_cdr.set_ylabel('Observed Freq', fontsize=9)
        ax_cdr.set_ylim(0, 1.18)
        ax_cdr.set_xticks(range(len(sorted_cdr)))
        ax_cdr.set_xticklabels(cdr_labels, fontsize=7.5, rotation=0)

        ax_fw.legend(bbox_to_anchor=(1.01, 0.5), loc='center left', ncol=2, fontsize=7.5)

        plt.tight_layout()
        plot_path = os.path.join(output_dir, f"{chain_name}_{gene_name}_Split_PWM.png")
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close(fig)

        image_paths[gene_name] = plot_path

    return image_paths
