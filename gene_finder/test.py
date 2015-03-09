# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 11:24:42 2014

@author: Meghan Tighe

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
    # TODO: implement this
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'
    else:
        return 'WTF?'

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
    # TODO: implement this
    result = ""
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
    >>> rest_of_ORF("ATGAATGTAG")
    'ATGAATGTAG'
    >>> rest_of_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCGAATG'
    """
    # TODO: implement this
    end_index = len(dna)
    i = 3
    stop_codons = ['TAG', 'TAA', 'TGA']
    while i <= len(dna)-2:
        if dna[i:i+3] in stop_codons:
            end_index = i
            break
        i+=3
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
    >>> find_all_ORFs_oneframe("ATGAATGTAG")
    ['ATGAATGTAG']
    >>> find_all_ORFs_oneframe("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG']
    """
    # TODO: implement this
    ORFs = []
    i = 0
    while i <len(dna):
        if dna[i:i+3] == "ATG":
            ORFs.append(rest_of_ORF(dna[i:]))
            i = len(rest_of_ORF(dna[i:])) + i
        else:
            i+=3
    return ORFs

def collapse(L):
    '''Converts a list of strings to a string by conncat-ing all elements of the list'''
    output = ""
    for s in L:
        output = output + s
    return output

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in all 3
        possible frames and returns them as a list.  By non-nested we mean that if an
        ORF occurs entirely within another ORF and they are both in the same frame,
        it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    >>> find_all_ORFs("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG']
    >>> find_all_ORFs("")
    []
    """
    # TODO: implement this
    ORFs = []
    for i in range(0,3):
        orfs_oneframe = find_all_ORFs_oneframe(dna[i:])
        ORFs.extend(orfs_oneframe)
    return ORFs


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    >>> find_all_ORFs_both_strands("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG', 'ATGCAT']
    """
    # TODO: implement this
    ORFs = []
    dna_rev_comp = get_reverse_complement(dna)
    ORFs.extend(find_all_ORFs(dna))
    ORFs.extend(find_all_ORFs(dna_rev_comp))
    return ORFs



def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    >>> longest_ORF("ATGCATGAATGTAG")
    'ATGCATGAATGTAG'
    """
    all_ORFs = find_all_ORFs_both_strands(dna)
    longest = max(all_ORFs, key=len)
    return longest
    # TODO: implement this


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    longest = ''
    i = 1
    while i<=num_trials:
        shuffled = shuffle_string(dna)
        a_long = longest_ORF(shuffled)
        if len(a_long)>len(longest):
            longest = a_long
        i+=1
    return len(longest)
    # TODO: implement this


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
    # TODO: implement this
    amino_acids= ''
    i=0
    while i<=len(dna)-3:
        amino_acids = amino_acids + (aa_table[dna[i:i+3]])
        i = i + 3
    return amino_acids


def gene_finder(dna):
    """ Returns the amino acid sequences coded by all genes that have an ORF
        larger than the specified threshold.
        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    threshold = longest_ORF_noncoding(dna, 1500)
    print "threshold is ", threshold
    orfs = find_all_ORFs_both_strands(dna)
    good_orfs_aas=[]
    for item in orfs:
        if len(item)>=threshold:
            good_orfs_aas.append(coding_strand_to_AA(item))
    return good_orfs_aas

def salmonella_gene_finder():
    dna = load_seq("./data/X73525.fa")
    genes = gene_finder(dna)
    for gene in genes:
        print gene


if __name__ == "__main__":
    import doctest
    doctest.testmod()
    salmonella_gene_finder()