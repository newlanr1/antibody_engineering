# src/reverse_complement.py

def reverse_string(seq):
    """Returns the reversed sequence using slicing."""
    return seq[::-1]

def complement(seq):
    """Returns the complementary sequence string using translation tables."""
    # str.maketrans is the modern, optimized Python 3 way to do this
    transtable = str.maketrans('acgtnACGTN', 'tgcanTGCAN')
    return seq.translate(transtable)

def reverse_complement(seq):
    """Returns the complete reverse complement of a DNA string."""
    reversed_seq = reverse_string(seq)
    return complement(reversed_seq)

if __name__ == "__main__":
    # This block allows you to test the script quickly 
    # without running it when imported elsewhere.
    sample_dna = "ATGCGTAN"
    
    print(f"Original Sequence:  {sample_dna}")
    print(f"Reversed:           {reverse_string(sample_dna)}")
    print(f"Complemented:       {complement(sample_dna)}")
    print(f"Reverse Complement: {reverse_complement(sample_dna)}")