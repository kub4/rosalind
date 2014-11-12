#!/usr/bin/env python3
"""
Introduction to the Bioinformatics Armory
=========================================

Given: A DNA string s of length at most 1000 bp.

Return: Four integers (separated by spaces) representing the respective
number of times that the symbols 'A', 'C', 'G', and 'T' occur in s. 

Sample Dataset
--------------

AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC

Sample Output
-------------

20 12 17 21

"""
import sys
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

# extract dna from the file, use the biopython sequence object
# and define its alphabet just for the joy of it
with open(sys.argv[1], 'r') as in_file:
  dna = Seq(''.join(in_file.read().upper().split()), generic_dna)

# define the order of bases in the output and create a counter list
baseorder  = "ACGT"
basecount = []

# count the bases
for base in baseorder:
  basecount.append(dna.count(base))

# print the results
print(" ".join([str(x) for x in basecount]))
