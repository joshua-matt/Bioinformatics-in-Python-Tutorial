from DNA.structures import *
from DNA.utilities import *
import random

# Procedure: validate_seq
# Purpose: Determine if a string represents a valid DNA sequence
# Parameters: dna_seq, a string
# Produces: uppercase of the string if only C, G, A, T characters (case insensitive), raises Exception otherwise
def validate_seq(dna_seq):
    """Determine if a string represents a valid DNA sequence"""
    seq_upper = dna_seq.upper()
    for nuc in seq_upper:
        if nuc not in NUCLEOTIDES:
            raise Exception(f"Invalid DNA string {dna_seq}")
    return seq_upper

# Procedure: random_seq
# Purpose: Generate a random DNA sequence
# Parameters: length, a positive integer
# Produces: a random string of C, G, T, A of length 'length'
def random_seq(length):
    """Generate a random DNA sequence"""
    return ''.join([random.choice(NUCLEOTIDES) for i in range(length)])

# Procedure: count_nuc
# Purpose: Count the number of each nucleotide in a DNA sequence
# Parameters: dna_seq, a string
# Produces: a dictionary with keys "C", "G", "T", "A" and values of their respective counts in dna_seq
def count_nuc(dna_seq):
    """Count the number of each nucleotide in a DNA sequence"""
    nuc_counts = {nuc: 0 for nuc in NUCLEOTIDES}
    for nuc in validate_seq(dna_seq):
        nuc_counts[nuc] += 1
    return nuc_counts

# Procedure: transcription
# Purpose: Transcribe a DNA string into the corresponding RNA string
# Parameters: dna_seq, a string
# Produces: rna_seq, dna_seq with all "T" replaced with "U"
def transcription(dna_seq):
    """Transcribe a DNA string into the corresponding RNA string, replacing Thymine with Uracil"""
    return validate_seq(dna_seq).replace("T", "U")

# Procedure: reverse_complement
# Purpose: Create the complement of the DNA string then reverse it
# Parameters: dna_seq, a string
# Produces: rev_comp, the reverse of dna_seq where "A" and "T" are swapped and "G" and "C" are swapped
def reverse_complement(dna_seq):
    """Swap Adenine and Thymine, Guanine and Cytosine, and reverse the result"""
    return ''.join([NUCLEOTIDE_COMPLEMENTS[nuc] for nuc in validate_seq(dna_seq)])[::-1]

# Procedure: display_seq
# Purpose: Show all steps of DNA analysis and processing
# Parameters: dna_seq, a string
# Produces: None, prints sequence, sequence length, nucleotide counts, DNA -> RNA transcription to console
def display_seq(dna_seq):
    """Show all steps of DNA analysis and processing"""
    valid_seq = validate_seq(dna_seq)
    printn(f"Sequence: {valid_seq}")
    printn(f"1) Sequence Length: {len(valid_seq)}")
    printn(f"2) Nucleotide Frequency: {count_nuc(valid_seq)}")
    printn(f"3) DNA/RNA Transcription: {transcription(valid_seq)}")
    print(f"4) DNA String + Reverse Complement:")
    print(f"5' {valid_seq} 3'")
    print(f"   {'|' * len(valid_seq)}")
    print(f"3' {reverse_complement(valid_seq)} 5'")