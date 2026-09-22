# src/fasta_sequence_lengths.py
import argparse
from Bio import SeqIO

def get_sequence_lengths(file_path):
    """Returns a dictionary of {sequence_id: length} from a FASTA file."""
    return {record.id: len(record.seq) for record in SeqIO.parse(file_path, "fasta")}

def find_extremes(sequences_dict):
    """Returns the IDs of the longest and shortest sequences in a dictionary."""
    if not sequences_dict:
        return None, None
    longest_id = max(sequences_dict, key=sequences_dict.get)
    shortest_id = min(sequences_dict, key=sequences_dict.get)
    return longest_id, shortest_id

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate sequence lengths from a FASTA file.")
    parser.add_argument("file_path", help="Path to the FASTA file")
    args = parser.parse_args()

    lengths = get_sequence_lengths(args.file_path)
    if not lengths:
        print("No sequences found.")
    else:
        long_id, short_id = find_extremes(lengths)
        for seq_id, length in lengths.items():
            print(f"{seq_id}: {length}")
        print(f"\nLongest: {long_id} ({lengths[long_id]}bp)")
        print(f"Shortest: {short_id} ({lengths[short_id]}bp)")