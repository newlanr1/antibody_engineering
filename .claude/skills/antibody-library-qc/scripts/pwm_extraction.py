"""
Parse a library "design rules" Excel workbook into per-gene PWM dictionaries.

Expected sheet layout (one block per gene/design, blocks stacked vertically):
  - A cell in column index 1 (0-indexed) holding the design/gene name marks the
    start of a block.
  - The row immediately below the name row holds the full wild-type amino acid
    sequence, one residue per cell, starting at column index 4.
  - Three rows below the name row ("WT row"), cells containing the literal
    string "WT" mark the start of a designed (e.g. CDR) position's mutation
    column. Each subsequent column to the right of a "WT" marker is another
    designed position, until a blank cell is hit.
  - For each designed position column, rows start+5 through start+24 (20 rows)
    list amino acid letters in column `col_idx` and their design frequency in
    column `col_idx + offset`. A frequency > 0 means that amino acid is an
    allowed substitution at that position.

This layout matches AffMat_design.xlsx. If your workbook uses a different
convention, adjust the row/column offsets below rather than the call sites.
"""

import pandas as pd


def extract_full_pwm(filepath, sheet_name):
    """Return {design_name: {'wt_sequence': str, 'pwm': {pos: {aa: freq}},
    'allowed_mutations': {pos: [aa, ...]}}} for every gene block in the sheet."""
    df = pd.read_excel(filepath, sheet_name=sheet_name)
    all_indices = df[df.iloc[:, 1].notna()].index.tolist()
    design_start_indices = all_indices[0::2]

    gene_pwms = {}

    for i in range(len(design_start_indices)):
        start_idx = design_start_indices[i]
        end_idx = design_start_indices[i + 1] if i + 1 < len(design_start_indices) else len(df)

        design_name = df.iloc[start_idx, 1]
        wt_aa_sequence = "".join(df.iloc[start_idx + 1, 4:].dropna().astype(str))

        wt_row_idx = start_idx + 3
        if wt_row_idx >= end_idx:
            continue

        wt_row = df.iloc[wt_row_idx]
        wt_marker_cols = wt_row[wt_row == 'WT'].index.tolist()

        pwm_dict = {}
        allowed_mutations = {}

        for col_name in wt_marker_cols:
            col_idx = df.columns.get_loc(col_name)
            offset = 1
            while (col_idx + offset) < len(df.columns):
                if pd.isna(df.iloc[wt_row_idx, col_idx + offset]):
                    break

                seq_position = (col_idx + offset) - 4

                aa_probs = {}
                allowed_aas = []
                for aa_row_offset in range(5, 25):
                    aa_label = str(df.iloc[start_idx + aa_row_offset, col_idx]).strip()
                    prob = df.iloc[start_idx + aa_row_offset, col_idx + offset]

                    if pd.notna(prob) and isinstance(prob, (int, float)):
                        aa_probs[aa_label] = float(prob)
                        if prob > 0:
                            allowed_aas.append(aa_label)

                if aa_probs:
                    pwm_dict[seq_position] = aa_probs
                if allowed_aas:
                    allowed_mutations[seq_position] = allowed_aas

                offset += 1

        gene_pwms[design_name] = {
            'wt_sequence': wt_aa_sequence,
            'pwm': pwm_dict,
            'allowed_mutations': allowed_mutations,
        }

    return gene_pwms
