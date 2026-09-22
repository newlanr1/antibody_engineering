# src/mutagenesis.py
import argparse
import csv
from src.codon_tables import (
    STANDARD_CODON_TABLE,
    PREFERRED_CODONS_YEAST,
    PREFERRED_CODONS_HUMAN,
    PREFERRED_CODONS_ECOLI
)

# Map string arguments to the correct dictionary
SPECIES_PREFERENCES = {
    'yeast': PREFERRED_CODONS_YEAST,
    'human': PREFERRED_CODONS_HUMAN,
    'ecoli': PREFERRED_CODONS_ECOLI
}

def _get_aa_for_codon(codon):
    for aa, codons in STANDARD_CODON_TABLE.items():
        if codon in codons:
            return aa
    return None

def generate_single_mutations(dna_seq, species='yeast'):
    """Generates single amino acid substitutions optimized for the given species."""
    preferred_codons = SPECIES_PREFERENCES[species]
    mutated_seqs = []
    
    for i in range(0, len(dna_seq), 3):
        original_aa = _get_aa_for_codon(dna_seq[i:i+3])
        if not original_aa: continue

        for aa, preferred_codon in preferred_codons.items():
            if aa != 'C' or original_aa == 'C':
                if aa != original_aa:
                    mutated = dna_seq[:i] + preferred_codon + dna_seq[i+3:]
                    mutated_seqs.append(mutated)
    return mutated_seqs

def generate_double_mutations(dna_seq, species='yeast'):
    """Generates pairwise double amino acid substitutions optimized for the given species."""
    preferred_codons = SPECIES_PREFERENCES[species]
    mutated_seqs = []
    length = len(dna_seq)
    
    for i in range(0, length, 3):
        for j in range(i + 3, length, 3):
            aa1 = _get_aa_for_codon(dna_seq[i:i+3])
            aa2 = _get_aa_for_codon(dna_seq[j:j+3])
            if not aa1 or not aa2: continue

            for new_aa1, pref1 in preferred_codons.items():
                if new_aa1 == 'C' and aa1 != 'C': continue
                for new_aa2, pref2 in preferred_codons.items():
                    if new_aa2 == 'C' and aa2 != 'C': continue
                    if new_aa1 != aa1 or new_aa2 != aa2: 
                        mutated = dna_seq[:i] + pref1 + dna_seq[i+3:j] + pref2 + dna_seq[j+3:]
                        mutated_seqs.append(mutated)
    return mutated_seqs

def generate_triple_mutations(dna_seq, species='yeast'):
    """Generates triple amino acid substitutions optimized for the given species."""
    preferred_codons = SPECIES_PREFERENCES[species]
    mutated_seqs = []
    length = len(dna_seq)
    
    for i in range(0, length, 3):
        for j in range(i + 3, length, 3):
            for k in range(j + 3, length, 3):
                aa1 = _get_aa_for_codon(dna_seq[i:i+3])
                aa2 = _get_aa_for_codon(dna_seq[j:j+3])
                aa3 = _get_aa_for_codon(dna_seq[k:k+3])
                if not aa1 or not aa2 or not aa3: continue

                for new_aa1, pref1 in preferred_codons.items():
                    if new_aa1 == 'C' and aa1 != 'C': continue
                    for new_aa2, pref2 in preferred_codons.items():
                        if new_aa2 == 'C' and aa2 != 'C': continue
                        for new_aa3, pref3 in preferred_codons.items():
                            if new_aa3 == 'C' and aa3 != 'C': continue
                            
                            # Ensure at least one mutation occurs across the three positions
                            if new_aa1 != aa1 or new_aa2 != aa2 or new_aa3 != aa3: 
                                mutated = (dna_seq[:i] + pref1 + 
                                           dna_seq[i+3:j] + pref2 + 
                                           dna_seq[j+3:k] + pref3 + 
                                           dna_seq[k+3:])
                                mutated_seqs.append(mutated)
    return mutated_seqs

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate mutagenesis libraries.")
    parser.add_argument("type", choices=["single", "double", "triple"], help="Type of mutation library")
    parser.add_argument("sequence", help="Input DNA sequence")
    parser.add_argument("output_csv", help="Path to save results")
    parser.add_argument("--species", choices=["yeast", "human", "ecoli"], default="yeast", help="Target species for codon optimization")
    args = parser.parse_args()

    seq_upper = args.sequence.strip().upper()
    
    if args.type == "single":
        results = generate_single_mutations(seq_upper, args.species)
    elif args.type == "double":
        results = generate_double_mutations(seq_upper, args.species)
    else:
        results = generate_triple_mutations(seq_upper, args.species)
        
    with open(args.output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Mutated Sequence'])
        for r in results:
            writer.writerow([r])
            
    print(f"Generated {len(results)} {args.type} mutations optimized for {args.species} and saved to {args.output_csv}")