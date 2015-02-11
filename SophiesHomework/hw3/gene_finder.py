# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 11:24:42 2014

@author: Sophie 
why is this so hard lol ;___;
"""

from amino_acids import aa, codons
from load import load_seq
from random import shuffle

dna = load_seq("./data/X73525.fa")
codonlookup = dict(zip(aa, codons)) #creates dictionary mapping AA to condons

def collapse(L):
    """ Converts a list of strings to a string by concatenating all elements of the list """
    
    output = ""
    for s in L:
        output = output + s
    return output

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).
        
        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment
    """
    index = 0
    res = []
    while index < len(dna):
        codon = dna[index:index+3]
        for k in codonlookup:
            if codon in codonlookup[k]:
                res.append(k)
        index += 3
    return collapse(res) #collapse at end of most functions to make strings, not lists
        
def coding_strand_to_AA_unit_tests (my_input, expected_output):
    """ Unit tests for the coding_strand_to_AA function """
    print 'Translates a codon to an amino acid'
    print 'input:', my_input, "expected output:", expected_output,"actual output:", coding_strand_to_AA(my_input)
    print '    '

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence
    
        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    """
    index = 0
    res = []
    while index < len(dna):
        N = dna[index]
        
        if N == 'A':
            res.append('T')            
        elif N == 'T':
            res.append('A')
        elif N == 'G':
            res.append('C')
        elif N == 'C':
            res.append('G')
        else:
            print "If you see me you messed something up really badly. ):"
        index += 1
    return collapse(res[::-1]) #Reverses this nucleotide chain

def get_reverse_complement_unit_tests(my_input, expected_output):
    """ Unit tests for the get_complement function """
    print 'Finds complementary nucleotides'
    print 'input:', my_input, "expected output:", expected_output,"actual output:", get_reverse_complement(my_input)
    print '   '

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start codon and returns
        the sequence up to but not including the first in frame stop codon.  If there
        is no in frame stop codon, returns the whole string.
        
        dna: a DNA sequence
        returns: the open reading frame represented as a string
    """
    index = 0
    while index < len(dna):
        codon = dna[index:index+3]
        index += 3
        if codon in ['TAG','TAA','TGA']:
            return dna[0:index-3]
    return dna

def rest_of_ORF_unit_tests(my_input ,expected_output): 
    """ Unit tests for the rest_of_ORF function """
    print 'Returns snippet of DNA up to, but not including stop codon'
    print 'input:', my_input, "expected output:", expected_output,"actual output:", rest_of_ORF(my_input)
    print '    '
    
def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence and returns
        them as a list.  This function should only find ORFs that are in the default
        frame of the sequence (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
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
    
  
def find_all_ORFs_oneframe_unit_tests(my_input, expected_output):
    """ Unit tests for the find_all_ORFs_oneframe function """
    print 'Finds all ORFs within one reading frame'
    print 'input:', my_input, "expected output:", expected_output,"actual output:", find_all_ORFs_oneframe(my_input)
    print '    '

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in all 3
        possible frames and returns them as a list.  By non-nested we mean that if an
        ORF occurs entirely within another ORF and they are both in                                             the same frame,
        it should not be included in the returned list of ORFs.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    """
    res = []
    p = 0
    while p < 3:
        dna_snippet = dna[p:len(dna)]
        res.append(collapse(find_all_ORFs_oneframe(dna_snippet)))
        p += 1
    return res

def find_all_ORFs_unit_tests(my_input, expected_output):
    """ Unit tests for the find_all_ORFs function """
    print 'Finds all ORFs within any of the 3 possible reading frames'
    print 'input:', my_input, "expected output:", expected_output,"actual output:", find_all_ORFs(my_input)
    print '    '
      
def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    """
    res = []
    dna_inverse = get_reverse_complement(dna)
    res.append(collapse(find_all_ORFs(dna)))
    res.append(collapse(find_all_ORFs(dna_inverse)))
    return res

def find_all_ORFs_both_strands_unit_tests(my_input, expected_output):
    """ Unit tests for the find_all_ORFs_both_strands function """
    print 'Finds all ORFs in both strands of DNA'
    print 'input:', my_input, "expected output:", expected_output,"actual output:", find_all_ORFs_both_strands(my_input)
    print '    '
 
def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string"""
    return max(find_all_ORFs_both_strands(dna),key=len)

def longest_ORF_unit_tests(my_input, expected_output):
    """ Unit tests for the longest_ORF function """
    print 'Returns longest ORF in both strands'
    print 'input:', my_input, "expected output:", expected_output,"actual output:", longest_ORF(my_input)
    print '    '
    
def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
#        of the specfied DNA sequence
#        
#        dna: a DNA sequence
#        num_trials: the number of random shuffles
#        returns: the maximum length longest ORF """
    i = 0
    res = []
    while i < num_trials:
        dna_shuffle = list(dna)
        shuffle(dna_shuffle)
        dna_shuffle = ''.join(dna_shuffle)
        res.append(len(longest_ORF(dna_shuffle))/9)
        i += 1
    return max(res)

print longest_ORF_noncoding(dna,1500)

def gene_finder(dna, threshold):
    """ Returns the amino acid sequences coded by all genes that have an ORF
        larger than the specified threshold.
        
        dna: a DNA sequence
        threshold: the minimum length of the ORF for it to be considered a valid
                   gene.
        returns: a list of all amino acid sequences whose ORFs meet the minimum
                 length specified.
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
     
print gene_finder(dna,720)

#Debugging unit tests
#Format is "Input", "Expected Output"
if __name__=="__main__":
    coding_strand_to_AA_unit_tests('ATGCCCGCTTT','MPA')
    get_reverse_complement_unit_tests('ATGCCCGCTTT','AAAGCGGGCAT')
    rest_of_ORF_unit_tests('ATGTGAA','ATG')
    find_all_ORFs_oneframe_unit_tests('ATGCATGAATGTAGATAGATGTGCCC','ATGCATGAATGTAGA, ATGTGCCC')
    find_all_ORFs_unit_tests('ATGCATGAATGTAG','ATGCATGAATGTAG, ATGAATGTAG, ATG')
    find_all_ORFs_both_strands_unit_tests('ATGCGAATGTAGCATCAAA', 'ATGCGAATG, ATGCTACATTCGCAT')
    longest_ORF_unit_tests('ATGCGAATGTAGCATCAAA','ATGCTACATTCGCAT')