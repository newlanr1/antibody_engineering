# src/dna_has_stop.py
import argparse

def has_stop_codon(dna_sequence, frame=0):
    """Checks if a DNA sequence has in-frame stop codons (TGA, TAG, TAA)."""
    stop_codons = {'tga', 'tag', 'taa'}
    dna_lower = dna_sequence.lower()
    
    for i in range(frame, len(dna_lower) - 2, 3):
        if dna_lower[i:i+3] in stop_codons:
            return True
    return False

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check for in-frame stop codons.")
    parser.add_argument("sequence", help="DNA sequence to check")
    parser.add_argument("--frame", type=int, default=0, choices=[0, 1, 2], help="Reading frame (0, 1, or 2)")
    args = parser.parse_args()

    if has_stop_codon(args.sequence, args.frame):
        print("Input sequence has an in-frame stop codon.")
    else:
        print("Input sequence has no in-frame stop codons.")