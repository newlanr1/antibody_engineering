"""
Two complementary ML/stats diagnostics for a library that's already been
scored by aligner.score_read_against_designs (i.e. has *_Design / *_True_Errors
/ *_Intended_Mutations columns, as produced by qc_pipeline.py):

- KL divergence between each designed position's *intended* amino acid
  distribution (from the PWM) and what was *actually observed* in sequencing.
  High divergence at a position means the library's real composition drifted
  from the design — e.g. oligo synthesis bias overrepresenting certain codons.

- ESM-2 pseudo-perplexity: a protein language model's view of how "natural"/
  stable a full sequence looks. Used here to check whether reads with more
  true framework errors also look less fit to the model — a sign that
  framework errors aren't just noise but are actually destabilizing.
"""

import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import entropy

AA_LIST = sorted(list("ACDEFGHIKLMNPQRSTVWY"))


def calculate_pwm_kl_divergence(df, rules_dict, chain_col, design_col, read_count_col, chain_name,
                                 best_window_col=None):
    """Return a DataFrame with one row per (gene, designed position): the KL
    divergence of the observed amino-acid distribution from the design's
    intended PWM at that position. Only positions present in a gene's `pwm`
    (i.e. designed/CDR positions, not framework) are scored — KL divergence
    against an intended distribution is meaningless for framework positions,
    which don't have one.

    If `best_window_col` is given (the window column aligner.score_reads_parallel
    already computed, e.g. "Heavy_Best_Window"), reuses that cached alignment
    instead of re-running the sliding-window search here. Omit it to fall
    back to recomputing from `chain_col` (e.g. standalone use, without
    qc_pipeline.py having run first)."""
    kl_results = []

    for gene_name, rules in rules_dict.items():
        pwm_target = rules['pwm']
        wt_seq = rules['wt_sequence']
        design_len = len(wt_seq)

        gene_reads = df[df[design_col] == gene_name]
        if len(gene_reads) == 0 or not pwm_target:
            continue

        pos_counts = {p: {aa: 0 for aa in AA_LIST} for p in pwm_target.keys()}
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
                for pos in pwm_target.keys():
                    if pos < len(best_window):
                        res = best_window[pos]
                        if res in pos_counts[pos]:
                            pos_counts[pos][res] += read_count

        if total_depth == 0:
            continue

        epsilon = 1e-6
        for pos, intended_dict in pwm_target.items():
            p_intended = np.array([intended_dict.get(aa, 0.0) for aa in AA_LIST])
            p_observed = np.array([pos_counts[pos][aa] for aa in AA_LIST]) / total_depth

            p_intended = np.clip(p_intended, epsilon, 1.0)
            p_intended /= p_intended.sum()

            p_observed = np.clip(p_observed, epsilon, 1.0)
            p_observed /= p_observed.sum()

            kl_div = entropy(p_observed, p_intended)

            kl_results.append({
                'Chain': chain_name,
                'Gene': gene_name,
                'Position': pos + 1,
                'WT_AA': wt_seq[pos] if pos < len(wt_seq) else '-',
                'KL_Divergence': round(float(kl_div), 4),
            })

    return pd.DataFrame(kl_results)


def load_esm2(model_name="facebook/esm2_t6_8M_UR50D"):
    """Load an ESM-2 tokenizer + masked-LM model once, to be reused across
    many compute_esm2_perplexity calls. Requires `torch` and `transformers`,
    and will download the model from Hugging Face on first use."""
    import torch
    from transformers import AutoTokenizer, EsmForMaskedLM

    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = EsmForMaskedLM.from_pretrained(model_name)
    model.eval()
    return tokenizer, model


def compute_esm2_perplexity(sequence, tokenizer, model):
    """Average sequence pseudo-perplexity under ESM-2. Lower = the model
    finds the sequence more natural/stable; higher = more surprising, which
    often tracks with destabilizing mutations or framework errors."""
    import torch

    if not isinstance(sequence, str) or len(sequence) == 0:
        return float("nan")

    inputs = tokenizer(sequence, return_tensors="pt")
    input_ids = inputs["input_ids"]

    with torch.no_grad():
        outputs = model(input_ids=input_ids)
        logits = outputs.logits

    log_probs = torch.log_softmax(logits, dim=-1)
    target_log_probs = log_probs[0, range(1, input_ids.shape[1] - 1), input_ids[0, 1:-1]]

    nll = -torch.mean(target_log_probs).item()
    return round(math.exp(nll), 3)


def score_esm2_perplexity_batch(sequences, model_name="facebook/esm2_t6_8M_UR50D"):
    """Score each unique sequence in `sequences` once (reads with identical
    sequences get identical scores, so dedup first to avoid redundant model
    calls) and return {sequence: perplexity}."""
    tokenizer, model = load_esm2(model_name)
    unique_seqs = pd.Series(sequences).dropna().unique()
    return {seq: compute_esm2_perplexity(seq, tokenizer, model) for seq in unique_seqs}


def compute_wt_reference_perplexity(rules_dict, model_name="facebook/esm2_t6_8M_UR50D"):
    """ESM-2 has no universal "this folds properly" perplexity cutoff — it's a
    relative, model-specific score, not a validated biophysical threshold. The
    one biologically grounded anchor available here is each design's own
    wild-type sequence: it's a real, working antibody scaffold by definition,
    so its perplexity is a reasonable reference for "at least as fold-compatible
    as the parental scaffold." Scores every design's wt_sequence and returns
    {'mean', 'min', 'max'} across designs for this chain — a chain usually has
    several designs with meaningfully different intrinsic perplexities (a read
    that's a *perfect* match to its own design's WT can still sit above the
    mean of other designs' WT scores, which isn't a bug — it's why the min/max
    range is reported too, not just a single line)."""
    tokenizer, model = load_esm2(model_name)
    wt_scores = [compute_esm2_perplexity(rules['wt_sequence'], tokenizer, model)
                 for rules in rules_dict.values()]
    if not wt_scores:
        return None
    return {"mean": sum(wt_scores) / len(wt_scores), "min": min(wt_scores), "max": max(wt_scores)}


# --- AbLang2: antibody-specific fitness LM (default backend) ---
#
# AbLang v1 (separate heavy/light models) is the architecture that would map
# most directly onto this pipeline's existing per-chain structure, but its
# weights are hosted at opig.stats.ox.ac.uk, which returns 403 Forbidden from
# this environment regardless of User-Agent/Referer -- not a transient issue,
# so it's not usable here. AbLang2's default checkpoint is paired-native
# (scores heavy+light jointly) and is hosted on Zenodo instead, which *is*
# reachable. To preserve the existing per-chain reporting shape, the same
# model is used to score each chain alone (passing "" for the other chain)
# -- see score_ablang2_chain_batch -- with the joint paired score kept as an
# additional, separate "pairing compatibility" metric rather than a
# replacement for per-chain fitness.
#
# AbLang2 also exposes a built-in pseudo_log_likelihood that masks and
# re-scores every position individually -- the textbook-correct way to
# compute pseudo-perplexity for a bidirectional LM, but benchmarked at
# ~8.5s/sequence here, i.e. ~100-200x slower than the single-forward-pass
# approximation below. At that rate even a modest 10k-row sample would take
# roughly a day, which is incompatible with this skill's whole approach to
# scale (see qc_pipeline.py/aligner.py's "Performance at scale" notes in
# SKILL.md). compute_ablang2_perplexity instead does ONE unmasked forward
# pass and reads off the model's own probability for each true residue --
# the same fast-but-approximate methodology already used for ESM-2 above,
# so results from the two backends are at least methodologically comparable.

def load_ablang2(model_to_use="ablang2-paired"):
    """Load the AbLang2 paired antibody LM once, to be reused across many
    compute_ablang2_perplexity calls. Requires `ablang2` (pip install
    ablang2); downloads weights from Zenodo on first use."""
    import ablang2
    return ablang2.pretrained(model_to_use=model_to_use, random_init=False, ncpu=1, device="cpu")


def compute_ablang2_perplexity(heavy_seq, light_seq, ablang_model):
    """Single-forward-pass pseudo-perplexity for a heavy and/or light chain
    sequence. Pass "" for whichever chain you don't have (e.g. to score one
    chain alone) -- AbLang2 collapses an empty side out of its input format
    rather than erroring on it."""
    import torch

    heavy_seq = heavy_seq if isinstance(heavy_seq, str) else ""
    light_seq = light_seq if isinstance(light_seq, str) else ""
    if not heavy_seq and not light_seq:
        return float("nan")

    seq_str = f"<{heavy_seq}>|<{light_seq}>".replace("<>", "")
    tokens = ablang_model.tokenizer([seq_str], pad=True, w_extra_tkns=False, device=ablang_model.used_device)

    with torch.no_grad():
        logits = ablang_model.AbLang(tokens)

    special = torch.Tensor(ablang_model.tokenizer.all_special_tokens).to(ablang_model.used_device)
    mask = ~torch.isin(tokens, special)
    idxs = mask.nonzero()
    log_probs = torch.log_softmax(logits, dim=-1)
    observed_logp = log_probs[0, idxs[:, 1], tokens[0, idxs[:, 1]]]
    nll = -observed_logp.mean().item()
    return round(math.exp(nll), 3)


def score_ablang2_chain_batch(sequences, ablang_model, chain_side):
    """Score each unique sequence in `sequences` alone (the other chain left
    blank). `chain_side` is 'heavy' or 'light' -- same shared model either
    way, just controls which slot the sequence goes into. Returns
    {sequence: perplexity}, mirroring score_esm2_perplexity_batch's shape."""
    unique_seqs = pd.Series(sequences).dropna().unique()
    scores = {}
    for seq in unique_seqs:
        if chain_side == "heavy":
            scores[seq] = compute_ablang2_perplexity(seq, "", ablang_model)
        else:
            scores[seq] = compute_ablang2_perplexity("", seq, ablang_model)
    return scores


def score_ablang2_paired_batch(heavy_sequences, light_sequences, ablang_model):
    """Joint heavy+light pairing-compatibility score: scores each unique
    (heavy, light) pair together as one combined sequence, rather than each
    chain alone. Returns {(heavy_seq, light_seq): perplexity}."""
    pairs = pd.DataFrame({"h": heavy_sequences, "l": light_sequences}).dropna().drop_duplicates()
    return {(h, l): compute_ablang2_perplexity(h, l, ablang_model)
            for h, l in zip(pairs["h"], pairs["l"])}


def compute_ablang2_wt_reference(rules_dict, ablang_model, chain_side):
    """AbLang2 analogue of compute_wt_reference_perplexity: {'mean', 'min', 'max'}
    perplexity of each design's own wild-type sequence, scored alone on
    `chain_side`. See compute_wt_reference_perplexity's docstring for why the
    range matters, not just the mean -- different designs have genuinely
    different intrinsic perplexities, so a perfect match to one design's own
    WT can legitimately sit above the mean of every other design's WT score."""
    scores = [
        compute_ablang2_perplexity(rules['wt_sequence'], "", ablang_model) if chain_side == "heavy"
        else compute_ablang2_perplexity("", rules['wt_sequence'], ablang_model)
        for rules in rules_dict.values()
    ]
    if not scores:
        return None
    return {"mean": sum(scores) / len(scores), "min": min(scores), "max": max(scores)}


def compute_ablang2_paired_wt_reference(heavy_rules_dict, light_rules_dict, ablang_model):
    """Joint-pairing WT reference: {'mean', 'min', 'max'} perplexity across every
    combination of a heavy design's WT sequence paired with a light design's WT
    sequence. Cheap -- tens to low hundreds of (gene, gene) combinations, not
    reads. The spread here is usually wider than either single-chain reference,
    since it compounds both chains' design-to-design variation."""
    scores = [
        compute_ablang2_perplexity(h_rules['wt_sequence'], l_rules['wt_sequence'], ablang_model)
        for h_rules in heavy_rules_dict.values()
        for l_rules in light_rules_dict.values()
    ]
    if not scores:
        return None
    return {"mean": sum(scores) / len(scores), "min": min(scores), "max": max(scores)}


def plot_perplexity_vs_errors(df, error_col, perplexity_col, intended_col, chain_name,
                               output_path, error_threshold=2, fold_reference=None,
                               model_label="Fitness LM"):
    """Scatter plot: true framework errors (x) vs fitness-LM perplexity (y), colored
    by intended-mutation count, with a vertical line at the QC error cutoff.
    Useful for spotting whether framework errors actually correlate with
    reduced predicted fitness, or whether they're fitness-neutral noise.

    `model_label` (e.g. "ESM-2", "AbLang2") is just the axis/title text —
    this function doesn't care which backend scored perplexity_col.

    `fold_reference`, if given, is a {'mean', 'min', 'max'} dict (see
    compute_wt_reference_perplexity / compute_ablang2_wt_reference /
    compute_ablang2_paired_wt_reference) describing the spread of perplexity
    across each design's own wild-type sequence. Different designs have
    genuinely different intrinsic perplexities, so this draws the full
    min-max range as a band (not just a single mean line) — a point sitting
    inside the band is as fold-compatible as *some* real design's own WT,
    even if it's above the mean of every other design's WT score, which is
    not a red flag. Only above the max is actually outside the WT range this
    library's own designs occupy."""
    plot_df = df.dropna(subset=[perplexity_col])

    fig, ax = plt.subplots(figsize=(8, 4.5))
    sns.scatterplot(
        data=plot_df, x=error_col, y=perplexity_col, hue=intended_col,
        palette='viridis', s=60, ax=ax,
    )
    ax.set_title(f'{chain_name}: True Errors vs. {model_label} Predicted Perplexity',
                 fontsize=11, fontweight='bold')
    ax.set_xlabel('True Framework Errors', fontsize=10)
    ax.set_ylabel(f'{model_label} Perplexity (Lower = Higher Fitness)', fontsize=10)
    ax.axvline(error_threshold, color='crimson', linestyle='--',
               label=f'QC Error Cutoff (≤{error_threshold})')

    if fold_reference is not None:
        ymin, ymax = ax.get_ylim()
        ax.axhspan(ymin, fold_reference["min"], color='forestgreen', alpha=0.18, zorder=0)
        ax.axhspan(fold_reference["min"], fold_reference["max"], color='forestgreen', alpha=0.08, zorder=0)
        ax.axhline(fold_reference["mean"], color='forestgreen', linestyle=':', linewidth=1.5,
                   label=f'Mean WT baseline ({fold_reference["mean"]:.3f})')
        ax.axhline(fold_reference["min"], color='forestgreen', linestyle='-', linewidth=0.8, alpha=0.6)
        ax.axhline(fold_reference["max"], color='forestgreen', linestyle='-', linewidth=0.8, alpha=0.6,
                   label=f'WT design range ({fold_reference["min"]:.3f}-{fold_reference["max"]:.3f})')
        ax.set_ylim(ymin, ymax)

    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', title='Intended Muts')
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return output_path
