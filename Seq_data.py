default_transcribe = {'A': 'A', 'T': 'U', 'G': 'G', 'C': 'C', 'N': 'N'}

default_reverse_transcribe = {'A': 'A', 'U': 'T', 'G': 'G', 'C': 'C', 'N': 'N'}

default_complement_dna = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

default_complement_rna = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

valid_dna = {'A', 'T', 'G', 'C', 'N'}

valid_rna = {'A', 'U', 'G', 'C', 'N'}

valid_protein = {
    'A', 'R', 'N', 'D', 'C', 'G', 'Q', 'E', 'H', 'I',
    'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'
}

default_translate = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
    'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

