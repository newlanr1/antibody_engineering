# src/compare_sequences.py
import argparse
import csv
import sys
from Bio import Align

def load_sequences_from_csv(filepath):
    """Reads a CSV (Col 1: Name, Col 2: Sequence) into a list of tuples."""
    sequences = []
    with open(filepath, 'r', newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) >= 2:
                sequences.append((row[0].strip(), row[1].strip()))
    return sequences

def configure_aligner():
    """Configures a PairwiseAligner to count absolute exact amino acid matches."""
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1.0
    aligner.mismatch_score = 0.0
    aligner.open_gap_score = 0.0
    aligner.extend_gap_score = 0.0
    return aligner

def calculate_percent_identity(query_seq, ref_seq, aligner):
    """Calculates percent identity between two sequences using a provided aligner."""
    if not query_seq or not ref_seq:
        return 0.0
    
    alignments = aligner.align(query_seq, ref_seq)
    matches = alignments[0].score
    max_len = max(len(query_seq), len(ref_seq))
    
    return (matches / max_len) * 100.0 if max_len > 0 else 0.0

def find_best_match(query_seq, references, aligner):
    """Finds the best reference match and percent identity for a single query."""
    best_ref_name = "None"
    best_pid = 0.0
    
    for r_name, r_seq in references:
        pid = calculate_percent_identity(query_seq, r_seq, aligner)
        if pid > best_pid:
            best_pid = pid
            best_ref_name = r_name
            
    return best_ref_name, best_pid

if __name__ == "__main__":
    # Command-line interface is preserved in the execution block
    parser = argparse.ArgumentParser(description="Compare query amino acid sequences to reference sequences and calculate percent identity.")
    parser.add_argument("-q", "--query", required=True, help="Path to the query CSV file")
    parser.add_argument("-r", "--reference", required=True, help="Path to the reference CSV file")
    parser.add_argument("-o", "--output", required=True, help="Path for the output CSV file")
    args = parser.parse_args()

    try:
        references = load_sequences_from_csv(args.reference)
        queries = load_sequences_from_csv(args.query)
    except FileNotFoundError as e:
        print(f"Error loading files: {e}")
        sys.exit(1)

    if not references:
        print("Error: No reference sequences found.")
        sys.exit(1)

    aligner = configure_aligner()

    with open(args.output, 'w', newline='') as out_file:
        writer = csv.writer(out_file)
        writer.writerow(["Query_Name", "Query_Sequence", "Closest_Reference", "Percent_Identity"])
        
        for q_name, q_seq in queries:
            best_ref_name, best_pid = find_best_match(q_seq, references, aligner)
            writer.writerow([q_name, q_seq, best_ref_name, f"{best_pid:.2f}"])
            
    print(f"Success! Output written to {args.output}")