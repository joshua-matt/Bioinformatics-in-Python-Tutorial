from DNA.DNAToolkit import *
with open("rosalind_dna.txt") as f:
    dna_seq = f.read().strip()
    counts = count_nuc(dna_seq)
    print(' '.join([str(counts[nuc]) for nuc in NUCLEOTIDES]))
