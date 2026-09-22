# src/count_fasta_entries.py
import argparse
import sys

def count_fasta_entries(file_path):
    """Counts the number of entries in a FASTA file by counting '>' headers."""
    count = 0
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                count += 1
    return count

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count entries in a FASTA file.")
    parser.add_argument("file_path", help="Path to the FASTA file")
    args = parser.parse_args()

    try:
        num_entries = count_fasta_entries(args.file_path)
        print(f"Number of entries: {num_entries}")
    except FileNotFoundError:
        print(f"Error: File not found at {args.file_path}")
        sys.exit(1)