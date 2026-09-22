# src/codon_tables.py

# The Standard Genetic Code applies to Yeast, Human, and E. coli[cite: 7, 11]
STANDARD_CODON_TABLE = {
    'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'N': ['AAT', 'AAC'], 'D': ['GAT', 'GAC'], 'C': ['TGT', 'TGC'],
    'E': ['GAA', 'GAG'], 'Q': ['CAA', 'CAG'], 'G': ['GGT', 'GGC', 'GGA', 'GGG'],
    'H': ['CAT', 'CAC'], 'I': ['ATT', 'ATC', 'ATA'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'K': ['AAA', 'AAG'], 'M': ['ATG'], 'F': ['TTT', 'TTC'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'W': ['TGG'], 'Y': ['TAT', 'TAC'], 'V': ['GTT', 'GTC', 'GTA', 'GTG']
}

PREFERRED_CODONS_YEAST = {
    'A': 'GCT', 'R': 'CGT', 'N': 'AAT', 'D': 'GAT', 'C': 'TGT',
    'E': 'GAA', 'Q': 'CAA', 'G': 'GGT', 'H': 'CAT', 'I': 'ATT',
    'L': 'CTT', 'K': 'AAA', 'M': 'ATG', 'F': 'TTT', 'P': 'CCT',
    'S': 'TCT', 'T': 'ACT', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTT'
}

PREFERRED_CODONS_HUMAN = {
    'A': 'GCC', 'R': 'AGA', 'N': 'AAC', 'D': 'GAC', 'C': 'TGC',
    'E': 'GAG', 'Q': 'CAG', 'G': 'GGC', 'H': 'CAC', 'I': 'ATC',
    'L': 'CTG', 'K': 'AAG', 'M': 'ATG', 'F': 'TTC', 'P': 'CCC',
    'S': 'AGC', 'T': 'ACC', 'W': 'TGG', 'Y': 'TAC', 'V': 'GTG'
}

PREFERRED_CODONS_ECOLI = {
    'A': 'GCG', 'R': 'CGT', 'N': 'AAC', 'D': 'GAT', 'C': 'TGC',
    'E': 'GAA', 'Q': 'CAG', 'G': 'GGC', 'H': 'CAC', 'I': 'ATC',
    'L': 'CTG', 'K': 'AAA', 'M': 'ATG', 'F': 'TTC', 'P': 'CCG',
    'S': 'AGC', 'T': 'ACC', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTG'
}