# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 11:24:42 2014

@author: SOPHIE :D

"""

# you may find it useful to import these variables (although you are not required to use them)
from amino_acids import aa, codons, aa_table
import random
from load import load_seq

codonlookup = aa_table


def shuffle_string(s):
    """ Shuffles the characters in the input string
        NOTE: this is a helper function, you do not have to modify this in any way """
    return ''.join(random.sample(s,len(s)))


def collapse(L):
    """ Converts a list of strings to a string by concatenating all elements of the list """

    output = ""
    for s in L:
        output = output + s
    return output


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    N = nucleotide
    if N == 'A':
        return('T')
    elif N == 'T':
        return('A')
    elif N == 'G':
        return('C')
    elif N == 'C':
        return('G')
    else:
        return "If you see me you messed something up really badly. ):"

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
    index = 0
    res = []
    while index < len(dna):
        N = dna[index]
        res.append(get_complement(N))
        index += 1
    return collapse(res[::-1]) #Reverses this nucleotide chain


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
    """
    # TODO: implement this
    index = 0
    last_index = 0
    while index < len(dna):
        codon = dna[index:index+3]
        index += 3
        if codon in ['TAG','TAA','TGA']:
            return dna[0:index-3]
    return dna


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
    index = 0
    res = []
    while index < len(dna):
        codon=dna[index:index + 3]
        if codon == "ATG":
            res.append(rest_of_ORF(dna[index:]))
            index += len(res[-1]) #Lists are cyclic, this returns last value of list
        else:
            index += 3
    return res

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
    res = []
    p = 0
    while p < 3:
        dna_snippet = dna[p:len(dna)]
        res.append(collapse(find_all_ORFs_oneframe(dna_snippet)))
        p += 1
    return res

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    res = []
    dna_inverse = get_reverse_complement(dna)
    res.append(collapse(find_all_ORFs(dna)))
    res.append(collapse(find_all_ORFs(dna_inverse)))
    return res


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    return max(find_all_ORFs_both_strands(dna),key=len)


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    i = 0
    res = []
    while i < num_trials:
        dna_shuffle = list(dna)
        shuffle(dna_shuffle)
        dna_shuffle = ''.join(dna_shuffle)
        res.append(len(longest_ORF(dna_shuffle))/9)
        i += 1
    return max(res)

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
    index = 0
    res = []
    while index < len(dna):
        codon = dna[index:index+3]
        if codon in codonlookup:
            res.append(codonlookup[codon])
        index += 3
    return collapse(res) #collapse at end of most functions to make strings, not lists


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna
        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    ORFS = find_all_ORFs_both_strands(dna)
    i = 0
    res = []
    while i < len(ORFS):
        if len(ORFS[i]) > threshold:
            res.append(coding_strand_to_AA(ORFS[i]))
            i += 1
        else:
            i += 1
    return res

# print gene_finder(dna,720)

if __name__ == "__main__":
    import doctest
    doctest.testmod()