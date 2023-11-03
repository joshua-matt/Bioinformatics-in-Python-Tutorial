from DNA.DNAToolkit import *
from utilities import *
dna_seq = read_file("rosalind_revc.txt")
print(reverse_complement(dna_seq))