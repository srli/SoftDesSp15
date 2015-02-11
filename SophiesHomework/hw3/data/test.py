# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 11:09:31 2014

@author: sophie
"""
from amino_acids import aa, codons


#mylist = range(0,21)

#values = dict(zip(mylist, codons))
codonlookup = dict(zip(aa, codons))
#print codonlookup
def reverse_codonlookup(d , v):
    for k in d:
        if v in d[k]:
            return k
    raise ValueError

print reverse_codonlookup(codonlookup,'TGG')

#index = 0
#while index < len(dna):
#    codon = dna[index:index+3]
#    if codon == codons[1]:
#        print aa[1]
#    elif codon == codons[2]:
#        print aa[2]
#    else:
#        print 'J'
#    index = index + 3

#

#this finds the parts of the string that are codons
def find_codons(dna):
    index = 0
    res = []
    while index < len(dna):
        codon = dna[index:index+3]
        i = 0
        while i < len(codons):        
            if codon in codons[i]:
                res.append(codon)
                i = 25
            elif i < 21:
                i = i + 1
            else:
                return "null"
        index = index + 3
    print collapse(res)

def testcoding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).
        
        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment
    """
    index = 0
    while index < len(dna):
        codon = dna[index + 2]
        for k in codonlookup:
            if codon in codonlookup[k]:
                return k
        index = index + 3

print testcoding_strand_to_AA('TGGGAA')


#def histogram(s):
#    d = dict()
#    for c in s:
#        if c not in d:
#            d[c] = 1
#        else:
#            d[c] += 1
#    return d
#
#h = histogram('hello')
#print h
def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start codon and returns
        the sequence up to but not including the first in frame stop codon.  If there
        is no in frame stop codon, returns the whole string.
        
        dna: a DNA sequence
        returns: the open reading frame represented as a string
    """
    
    # YOUR IMPLEMENTATION HERE
    index = 0
    while index < len(dna):
        codon = dna[index:index+3]
        index = index + 3
        if codon in ['TAG','TAA','TGA']:
            return dna[0:index-3]
          
print rest_of_ORF('ATGTGAA')

     
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
            res.append(rest_of_ORF(dna))
            index = index + 3
        else:
            index = index + 3
    print res
    
find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")      
