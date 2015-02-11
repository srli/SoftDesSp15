# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 11:24:42 2014

@author: David Abrahams

"""

# you may find it useful to import these variables (although you are not required to use them)
from amino_acids import aa, codons, aa_table
import random
from load import load_seq

def shuffle_string(s):
    """ Shuffles the characters in the input string
        NOTE: this is a helper function, you do not have to modify this in any way """
    return ''.join(random.sample(s,len(s)))

### YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """

    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'
    elif nucleotide == 'T':
        return 'A'


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """

    result = ''

    for c in reversed(dna):
        result += get_complement(c)

    return result



def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start codon and returns
        the sequence up to but not including the first in frame stop codon.  If there
        is no in frame stop codon, returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGTGCCC")
    'ATGTGCCC'
    """
    i = 3
    end_index = len(dna)
    stop_codons = ['TAG', 'TAA', 'TGA']

    while i < len(dna) - 2:
        codon = dna[i:i+3]
        if codon in stop_codons:
            end_index = i
            break
        i += 3

    return dna[0:end_index]


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence and returns
        them as a list.  This function should only find ORFs that are in the default
        frame of the sequence (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    orfs = []
    start_codon = 'ATG'
    i = 0
    while i < len(dna) - 2:
        if dna[i:i+3] == start_codon:
            orf = rest_of_ORF(dna[i:])
            orfs.append(orf)
            i += len(orf)
        else:
            i += 3

    return orfs


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in all 3
        possible frames and returns them as a list.  By non-nested we mean that if an
        ORF occurs entirely within another ORF and they are both in the same frame,
        it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """

    orfs = [];

    for i in range(0, 3):
        orf_in_frame = find_all_ORFs_oneframe(dna[i:])
        orfs.extend(orf_in_frame)

    return orfs

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """

    orfs = []
    orfs.extend(find_all_ORFs(dna))
    orfs.extend(find_all_ORFs(get_reverse_complement(dna)))

    return orfs



def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    orfs = find_all_ORFs_both_strands(dna)
    longest = max(orfs, key=len)
    return longest


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    lengths = []
    for i in range(0, num_trials):
        shuffle = shuffle_string(dna)
        length = len(longest_ORF(shuffle))
        lengths.append(length)

    return max(lengths)


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """

    i = 0
    amino_acids = ''
    while i < len(dna) - 2:
        amino_acids += aa_table[dna[i:i+3]]
        i += 3
    return amino_acids



def gene_finder(dna):
    """ Returns the amino acid sequences coded by all genes that have an ORF
        larger than the specified threshold.

        dna: a DNA sequence
        threshold: the minimum length of the ORF for it to be considered a valid
                   gene.
        returns: a list of all amino acid sequences whose ORFs meet the minimum
                 length specified.
    """
    threshold = longest_ORF_noncoding(dna, 1500)
    all_orfs = find_all_ORFs_both_strands(dna)
    aa_sequences = []

    for orf in all_orfs:
        if (len(orf) >= threshold):
            aa_sequence = coding_strand_to_AA(orf)
            aa_sequences.append(aa_sequence)

    return aa_sequences

def salmonella_dna_finder():
    dna = load_seq("./data/X73525.fa")
    genes = gene_finder(dna)
    for gene in genes:
        print gene

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    salmonella_dna_finder()