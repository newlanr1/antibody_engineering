"""
Sliding-window aligner that scores an observed read against every wild-type
design scaffold and classifies each mismatch as an intended (designed) amino
acid substitution or a true error (unintended mismatch, e.g. synthesis or
cloning error, primer scarring, sequencing noise).

Sliding the full design window across the read (rather than assuming the read
starts exactly at the designed region) means untrimmed adapter/primer bases
or carried-over constant-domain sequence at either end don't get miscounted
as errors — the best-fit window absorbs them.

Each read is scored independently of every other read, so this is
embarrassingly parallel — score_reads_parallel fans the work out across CPU
cores. It's the dominant cost in the whole skill (it tests every read against
every design to find the best match), so downstream steps that already know
a read's winning design (plotting.py, bias_and_fitness.py) should reuse the
`best_window` this returns rather than re-deriving it from scratch.
"""

import os
from multiprocessing import Pool

import pandas as pd


def score_read_against_designs(read_seq, rules_dict):
    """Find the best-matching design for one read.

    `rules_dict` is {design_name: {'wt_sequence': str,
    'allowed_mutations': {pos: [aa, ...]}}} as produced by
    pwm_extraction.extract_full_pwm.

    Returns a 4-element pandas Series: [best_design_name, intended_mutation_count,
    true_error_count, best_window]. `best_window` is the actual substring of
    the read that produced the winning alignment — callers that already know
    which design a read belongs to can slice amino acids out of it directly
    instead of re-running the sliding-window search. Reads shorter than every
    design, or non-string input, get a sentinel result rather than raising.
    """
    if not isinstance(read_seq, str):
        return pd.Series(["No sequence", 0, None, None])

    best_design = "Unknown"
    min_errors = float('inf')
    associated_intended = 0
    best_window = None

    for design_name, rules in rules_dict.items():
        wt_seq = rules['wt_sequence']
        allowed_muts = rules['allowed_mutations']
        design_len = len(wt_seq)

        if len(read_seq) < design_len:
            continue

        for i in range(len(read_seq) - design_len + 1):
            window = read_seq[i:i + design_len]
            intended = 0
            errors = 0

            for pos in range(design_len):
                if window[pos] != wt_seq[pos]:
                    if pos in allowed_muts and window[pos] in allowed_muts[pos]:
                        intended += 1
                    else:
                        errors += 1

            if errors < min_errors:
                min_errors = errors
                associated_intended = intended
                best_design = design_name
                best_window = window

    if min_errors == float('inf'):
        return pd.Series(["Un-alignable", 0, None, None])

    return pd.Series([best_design, associated_intended, min_errors, best_window])


# --- Parallel scoring across CPU cores ---
#
# rules_dict is handed to each worker process once, via the Pool initializer,
# rather than re-pickled on every single read (which is what would happen if
# we instead did pool.map(lambda s: score_read_against_designs(s, rules_dict), ...)).

_WORKER_RULES = None


def _init_worker(rules_dict):
    global _WORKER_RULES
    _WORKER_RULES = rules_dict


def _score_one(read_seq):
    return score_read_against_designs(read_seq, _WORKER_RULES)


def score_reads_parallel(sequences, rules_dict, n_processes=None, chunksize=None):
    """Score every sequence in `sequences` against `rules_dict`, in parallel
    across CPU cores, preserving order. Returns a list of the same per-read
    results score_read_against_designs produces.

    Falls back to running in-process (no Pool) for small inputs or
    n_processes=1, since process startup isn't worth it below a few hundred
    reads. `chunksize` controls how many reads each worker grabs per task;
    left as None, it's sized automatically from the input length and core
    count to balance IPC overhead against load balancing.
    """
    n_processes = n_processes or os.cpu_count() or 1
    sequences = list(sequences)
    n = len(sequences)

    if n_processes <= 1 or n < 500:
        return [score_read_against_designs(s, rules_dict) for s in sequences]

    if chunksize is None:
        chunksize = max(1, n // (n_processes * 4))

    with Pool(processes=n_processes, initializer=_init_worker, initargs=(rules_dict,)) as pool:
        results = pool.map(_score_one, sequences, chunksize=chunksize)

    return results
