# src/count_sequence_repeats.py
import argparse
import sys
from collections import defaultdict
from Bio import SeqIO

def load_fasta_sequences(file_path):
    """Reads a FASTA file using Biopython and returns a dictionary of sequences."""
    sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def count_kmers(sequences, k):
    """Counts occurrences of k-mers (repeats of length k) across a dictionary of sequences."""
    repeat_counts = defaultdict(int)
    for sequence in sequences.values():
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            repeat_counts[kmer] += 1
    return repeat_counts

def find_most_frequent_repeat(repeat_counts):
    """Finds the most frequent repeat in the counts dictionary."""
    if not repeat_counts:
        return None, 0
    most_frequent = max(repeat_counts, key=repeat_counts.get)
    return most_frequent, repeat_counts[most_frequent]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count sequence repeats (k-mers) in a FASTA file.")
    parser.add_argument("file_path", help="Path to the FASTA file")
    parser.add_argument("n", type=int, help="Length of the sequence repeats to count")
    args = parser.parse_args()

    try:
        sequences = load_fasta_sequences(args.file_path)
    except FileNotFoundError:
        print(f"Error: File not found at {args.file_path}")
        sys.exit(1)

    if not sequences:
        print("No sequences found.")
        sys.exit(1)

    repeat_counts = count_kmers(sequences, args.n)
    
    if not repeat_counts:
        print(f"No repeats of length {args.n} found.")
        sys.exit(1)

    most_frequent, count = find_most_frequent_repeat(repeat_counts)
    
    print(f"Repeat counts of length {args.n}:")
    for repeat, cnt in repeat_counts.items():
        print(f"{repeat}: {cnt}")
    print(f"\nMost frequent repeat: {most_frequent} (Count: {count})")