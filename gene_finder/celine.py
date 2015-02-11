# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 11:24:42 2014

@author: YOUR NAME HERE

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
    >>> get_complement('a') # checks edge case behavior for catch-all case (else statement), note else statement should never be triggered
    '!'
    >>> get_complement('-')
    '!'
    """
    if nucleotide=='A':
        return 'T'
    elif nucleotide=='T':
        return 'A'
    elif nucleotide=='C':
        return 'G'
    elif nucleotide=='G':
        return 'C'
    else:
        return '!'

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'

    No other tests necessary, length and correct complement are already covered by the provided doctest
    """
    i=len(dna)-1
    rev_comp=""
    # alternately, for i in range (0,len(dna)): rev_comp+=get_complement(dna[len(dna)-i])
    while i >=0:
        rev_comp+=get_complement(dna[i])
        i-=1
    return rev_comp

def rest_of_terminated_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start codon and returns
        the sequence up to but not including the first in frame stop codon.  If there
        is no in frame stop codon, returns the whole string.

        NOTE: ANY PARTIAL CODONS WILL BE REMOVED FROM RETURNED STRING, rest_of_ORF will return strings as well that did not end with a stop codon

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    """
    dna_rest=""
    for i in range (0, (len(dna)/3)+1): # Note floor division for integer comparison, but shouldn't matter either way, +1 to get to last complete codon
    #alt: while i < len(dna)/3
        # cuts and reads next codon
        codon=dna[i*3:i*3+3] #index = i*3, since every third character starts the next codon
        # checks for stop codons TAG, TAA, TGA
        if codon=="TAG" or codon=="TAA" or codon=="TGA":
            return dna_rest
        else:
            dna_rest+=codon
    return dna_rest

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start codon and returns
        the sequence up to but not including the first in frame stop codon.  If there
        is no in frame stop codon, returns the whole string.

        NOTE: return statement includes strings that were not followed by a stop codon

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGAGACCC") # Edge case without stop codon
    'ATGAGACCC'
    """
    # orf = dna[0:3]
    # i = 3
    # while i < len(dna):
    #     if (dna[i:i+3] == 'TAG') or (dna[i:i+3] == 'TAA') or (dna[i:i+3] == 'TGA'):
    #         break
    #     orf += dna[i:i+3]
    #     i += 3
    # return orf

    dna_rest=""
    stop_index=len(dna)
    for i in range (len(dna)): # Note floor division for integer comparison, but shouldn't matter either way
        # cuts and reads next codon
        codon=dna[i:i+3] #index = i*3, since every third character starts the next codon
        # checks for stop codons TAG, TAA, TGA
        if codon=="TAG" or codon=="TAA" or codon=="TGA":
            stop_index=i
    return dna[0:stop_index]

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
    >>> find_all_ORFs_oneframe("TTATGGAATATG") #Test if it reads ORFs with start codon NOT on indices that are multiples of 3, but picks up later one
    ['ATG']
    """
    # Initialize list of ORFS
    l=[]
    dnaLeft=dna
    i=0
    # Find ATG start codon
    while i<len(dna):
        if dna[i:i+3]=='ATG':
            orf=rest_of_ORF(dnaLeft[i:])
            l.append(orf)
            i=i+len(orf)
        else:
            i+=3
    return l

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in all 3
        possible frames and returns them as a list.  By non-nested we mean that if an
        ORF occurs entirely within another ORF and they are both in the same frame,
        it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        NOTE: An alternative is to combine find_all_ORFs with find_all_ORFs_oneframe. Instead of splitting into three reading frames, find indices for all ATG start codons and find corresponding ORF

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    # loop through dna, with starting offsets of 0, 1, 2 (each offset corresponds to a different reading frame)
    orf_list=[]
    for offset in range (0,3): #range from 0-3 to get 0,1,2
        orf_list.extend(find_all_ORFs_oneframe(dna[offset:]))
    return orf_list

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    # get reverse complement
    complement=get_reverse_complement(dna)
    # extend list to include complement's ORFs
    orf_list=find_all_ORFs(dna)
    orf_list.extend(find_all_ORFs(complement))
    return orf_list


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    orfs=find_all_ORFs_both_strands(dna)
    long_orf=orfs[0]
    for frame in orfs:
        if len(frame)>len(long_orf):
            long_orf=frame
    return long_orf


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    record_orf=""
    for i in range(num_trials):
        shuffle=shuffle_string(dna)
        test_orf=longest_ORF(shuffle)
        if len(test_orf)>len(record_orf):
            record_orf=test_orf
    return len(record_orf)

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
    # initialize aa strand
    aa_strand=""
    # read codons
    for codon_index in range(0,len(dna),3):
        if len(dna)-codon_index>2:
            cod=dna[codon_index:codon_index+3]
            amino_acid=aa_table[cod]
            aa_strand+=amino_acid
    return aa_strand

def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    aa_list=[]
    # calculate threshold, arbitrarily set to 95% of len(longest_orf_noncoding)
    threshold=longest_ORF_noncoding(dna,1500)
    orfs=find_all_ORFs_both_strands(dna)
    for each in orfs:
        if len(each)>threshold:
            aa_list.append(coding_strand_to_AA(each))
    return aa_list

if __name__ == "__main__":
    #import doctest
    #doctest.testmod()
    from load import load_seq
    dna = load_seq("./data/X73525.fa")
    res= gene_finder(dna)
    print len(res)
    print res
d