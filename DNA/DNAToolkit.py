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
    # Pythonic approach, slightly faster
    # mapping = str.maketrans("ATCG, "TAGC")
    # return validate_seq(dna_seq).translate(mapping)[::-1]

# Procedure: gc_content
# Purpose: Compute percentage of sequence that is guanine or cytosine
# Parameters: dna_seq, a string
#             validate=True, whether to validate dna_seq
#             precision=6, number of decimals in GC content
# Produces: perc, a float in [0,100] corresponding to the percentage of characters in string that are 'G' or 'C'
def gc_content(dna_seq, validate=True, precision=6):
    """Compute percentage of sequence that is guanine or cytosine"""
    valid_seq = dna_seq
    if validate:
        valid_seq = validate_seq(dna_seq)
    return round(((valid_seq.count("G") + valid_seq.count("C"))/len(valid_seq)) * 100, precision)

# Procedure: gc_content_subseq
# Purpose: Compute GC contents of length k subsequences. k=20 by default
# Parameters: dna_seq, a string
#             k=20, the window size
# Produces: percs, a list of the GC contents of the subsequences
def gc_content_subseq(dna_seq, k=20):
    """Compute GC contents of length k subsequences of RNA/DNA. k=20 by default"""
    valid_seq = validate_seq(dna_seq)
    return [gc_content(valid_seq[i:i+k], validate=False) for i in range(0, len(dna_seq) - k + 1, k)]

# Procedure: translate_seq
# Purpose: Translate a DNA sequence into the amino acid sequence it codes for
# Parameters: dna_seq, a string
#             init_pos, where to start translating
# Produces: translated, a list of one-letter amino acid codes
def translate_seq(dna_seq, init_pos=0):
    """Translate a DNA sequence into the amino acid sequence it codes for"""
    valid_seq = validate_seq(dna_seq)
    return [DNA_CODONS[valid_seq[i:i+3]] for i in range(init_pos, len(valid_seq)-2, 3)]

# Procedure: codon_usage
# Purpose: Calculate the proportions of each codon encoding a particular amino acid
# Parameters: dna_seq, a string
#             amino, the one-letter code for an amino acid
# Produces: usage, a dictionary with amino acid keys and float values
def codon_usage(dna_seq, amino):
    """Calculate the proportions of each codon encoding a particular amino acid"""
    codons = []
    valid_seq = validate_seq(dna_seq)
    for pos in range(0, len(valid_seq)-2, 3):
        if DNA_CODONS[valid_seq[pos:pos+3]] == amino:
            codons.append(valid_seq[pos:pos+3])

    usage = dict(Counter(codons))
    total = len(codons)
    for seq in usage:
        usage[seq] = round(usage[seq] / total, 2)
    return usage

### My methods

# Procedure: random_seq
# Purpose: Generate a random DNA sequence
# Parameters: length, a positive integer
# Produces: a random string of C, G, T, A of length 'length'
def random_seq(length):
    """Generate a random DNA sequence"""
    return ''.join([random.choice(NUCLEOTIDES) for i in range(length)])

# Procedure: display_seq
# Purpose: Show all steps of DNA analysis and processing
# Parameters: dna_seq, a string
# Produces: None, prints sequence, sequence length, nucleotide counts, DNA -> RNA transcription to console
def display_seq(dna_seq):
    """Show all steps of DNA analysis and processing"""
    valid_seq = validate_seq(dna_seq)
    revc = reverse_complement(valid_seq)
    printn(f"Sequence: {valid_seq}")
    printn(f"1) Sequence Length: {len(valid_seq)}")
    printn(f"2) Nucleotide Frequency: {count_nuc(valid_seq)}")
    printn(f"3) DNA/RNA Transcription: {transcription(valid_seq)}")
    print(f"4) DNA String + Reverse Complement:")
    print(f"5' {valid_seq} 3'")
    print(f"   {'|' * len(valid_seq)}")
    print(f"3' {revc[::-1]} 5' (Complement)")
    printn(f"5' {revc} 3' (Reverse Complement)")
    printn(f"5) GC Content: {gc_content(valid_seq)}%")
    printn(f"6) GC Content in Subsequences k=5: {gc_content_subseq(valid_seq, k=5)}")
    printn(f"7) Amino Acid Sequence: {translate_seq(valid_seq)}")
    printn(f"8) Codon Usage (L): {codon_usage(valid_seq, 'L')}")
