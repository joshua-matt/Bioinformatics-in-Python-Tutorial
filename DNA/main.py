from DNA.DNAToolkit import *

dna_seq1 = "CAGTGGCA" # CAGTGGCA
dna_seq2 = "CAGTGGBA" # False
dna_seq3 = "" #
dna_seq4 = "cAgTGgcA" # CAGTGGCA

dna_seqs = [dna_seq1, dna_seq2, dna_seq3, dna_seq4]

for seq in dna_seqs:
    print(validate_seq(seq))

print(random_seq(5))

print(count_nuc(dna_seq3))

display_seq(random_seq(50))
