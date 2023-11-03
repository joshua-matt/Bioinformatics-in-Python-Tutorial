from DNA.DNAToolkit import *
from utilities import *
dna_seq = read_file("rosalind_rna.txt")
print(transcription(dna_seq))