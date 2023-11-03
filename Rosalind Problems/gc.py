from DNA.DNAToolkit import *
from utilities import *

def max_GC_FASTA(data):
    max_key = list(data.keys())[0]
    max_value = 0

    for key in data:
        gc = gc_content(data[key])
        if gc > max_value:
            max_key = key
            max_value = gc

    print(max_key)
    print(max_value)

max_GC_FASTA(read_FASTA("rosalind_gc.txt"))