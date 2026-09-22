# src/fastq_to_fasta.py
import argparse
from Bio import SeqIO

def convert_fastq_to_fasta(fastq_file, fasta_file):
    """Converts a FASTQ file to a FASTA file using Biopython."""
    count = SeqIO.convert(fastq_file, "fastq", fasta_file, "fasta")
    return count

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert FASTQ to FASTA.")
    parser.add_argument("input_fastq", help="Path to input .fastq file")
    parser.add_argument("output_fasta", help="Path to output .fasta file")
    args = parser.parse_args()
    
    records_converted = convert_fastq_to_fasta(args.input_fastq, args.output_fasta)
    print(f"Conversion complete: {records_converted} records converted to {args.output_fasta}")