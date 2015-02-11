# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 11:24:42 2014

@author: YOUR NAME HERE

"""

# you may find it useful to import these variables (although you are not required to use them)
from amino_acids import aa, codons, aa_table
import random
#importing regex stuff
import re
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
    >>> get_complement('U')
    ''
    >>> get_complement('Z')
    ''
    """
    # TODO: implement this
    if nucleotide=='A':
        return 'T'
    elif nucleotide=='T':
        return 'A'
    elif nucleotide=='G':
        return 'C'
    elif nucleotide=='C':
        return 'G'
    else:
        return ""

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
    res = "";
    for c in dna[::-1]:
        res = res + get_complement(c);
    return res;

def rest_of_ORF(dna):
    """
    I added the third unit test to see what would happen if there was no stop codon

    Takes a DNA sequence that is assumed to begin with a start codon and returns
        the sequence up to but not including the first in frame stop codon.  If there
        is no in frame stop codon, returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGTGCC")
    'ATGTGCC'
    >>> rest_of_ORF("ATGGCCTAG")
    'ATGGCC'
    """
    # TODO: implement this
    end_codons = ['TAG', 'TAA', 'TGA']
    current_codon = "ATG"
    index = 0;
    while current_codon not in end_codons and index < len(dna):
        index = index + 3;
        try:
            current_codon = dna[index:index+3]
        except IndexError:
            break;
    return dna[:index]

def find_all_ORFs_oneframe(dna):
    """
    I added a second unit test to experiment with multiple frames in one sequence

    Finds all non-nested open reading frames in the given DNA sequence and returns
        them as a list.  This function should only find ORFs that are in the default
        frame of the sequence (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGCATGAATAGATGGCCCGGTAGATGTGTTTTAAACCCGGGTAG")
    ['ATGCATGAA', 'ATGGCCCGG', 'ATGTGTTTTAAACCCGGG']
    """
    # TODO: implement this
    #print dna
    stop_codons = ['TAG', 'TGA', 'TAA']
    opened = False
    orf = ""
    orfs = []
    for i in range(0, len(dna)):
        if dna[i:i+3]=='ATG' and i%3==0 and opened==False:
            opened=True
        if opened==True:
            if dna[i:i+3] in stop_codons and i%3==0:
                opened=False
                orfs.append(orf)
                orf=""
            else:
                orf = orf + dna[i]
    if orf!="":
        orfs.append(orf)
        orf=""
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
    # TODO: implement this
    orfs = []
    for i in range(3):
        orfs.extend(find_all_ORFs_oneframe(dna[i:]))
    return orfs
    #start_indexes = [m.start() for m in re.finditer('ATG', dna)]
    #orfs = []
    #for index in start_indexes:
    #    orfs.append(rest_of_ORF(dna[index:]))
    #return orfs


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    >>> find_all_ORFs_both_strands("GGGAAACCCGGG")
    []
    >>> find_all_ORFs_both_strands("")
    []
    """
    # TODO: implement this
    strand2 = get_reverse_complement(dna)
    temp = find_all_ORFs(dna)
    temp.extend(find_all_ORFs(strand2))
    return temp


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    >>> longest_ORF("")
    ''
    >>> longest_ORF(2)
    ''
    """
    # TODO: implement this
    try:
        return max(find_all_ORFs_both_strands(dna), key=len);
    except ValueError:
        return ""
    except TypeError:
        return ""


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF

        Unit Test: As we run more tests, we get the longest orf possible using the given proteins

        >>> longest_ORF_noncoding('ATGCGAATGTAGCATCAAA', 1000)
        19
        >>> longest_ORF_noncoding('ATGCATGAATAGATGGCCCGGTAGATGTGTTTTAAACCCGGGTAG', 1000)
        45
        """
    # TODO: implement this
    strand = dna
    max_len = len(longest_ORF(strand))
    max_orf = ""
    for i in range(num_trials):
        strand = "".join(random.sample(strand, len(strand)))
        if len(longest_ORF(strand)) > max_len:
            max_len = len(longest_ORF(strand))
            max_orf = longest_ORF(strand)
    return max_len

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
    acids = ""
    codons = []
    for i in range(0, len(dna)-2, 3):
        codons.append(dna[i:i+3])
    acids = "".join(aa_table[i] for i in codons)
    return acids
    #pass

def gene_finder(dna, threshold):
    """ Returns the amino acid sequences coded by all genes that have an ORF
        larger than the specified threshold.

        dna: a DNA sequence
        threshold: the minimum length of the ORF for it to be considered a valid
                   gene.
        returns: a list of all amino acid sequences whose ORFs meet the minimum
                 length specified.
    >>> gene_finder('ATGCGAATGTAGCATCAAA', 3)
    ['MRM', 'MLHSH']
    >>> gene_finder('ATGCCCGCTTT', 3)
    ['MPA']
    >>> gene_finder('ATGCATGAATAGATGGCCCGGTAGATGTGTTTTAAACCCGGGTAG', 24)
    ['MNRWPGRCVLNPG']
    >>> gene_finder('ATGCATGAATAGATGGCCCGGTAGATGTGTTTTAAACCCGGGTAG', 12)
    ['MCFKPG', 'MNRWPGRCVLNPG']
    >>> gene_finder('ATGCATGAATAGATGGCCCGGTAGATGTGTTTTAAACCCGGGTAG', 4)
    ['MHE', 'MAR', 'MCFKPG', 'MNRWPGRCVLNPG', 'MH']
    """
    return [coding_strand_to_AA(i) for i in find_all_ORFs_both_strands(dna) if len(i) >= threshold]

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    data = load_seq("./data/X73525.fa")
    print gene_finder(data, longest_ORF_noncoding(data, 1500))