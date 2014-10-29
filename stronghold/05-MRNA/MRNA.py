#!/usr/bin/env python3

"""
Inferring mRNA from Protein
===========================

Given: A protein string of length at most 1000 aa.

Return: The total number of different RNA strings from which the protein could
have been translated, modulo 1,000,000. (Don't neglect the importance of the
stop codon in protein translation.)

Sample Dataset
--------------

MA

Sample Output
-------------

12

"""
import sys

with open(sys.argv[1], 'r') as in_file:
  protein = ''.join(in_file.read().split())

code_reversed = {
"A" : 4,
"C" : 2,
"D" : 2,
"E" : 2,
"F" : 2,
"G" : 4,
"H" : 2,
"I" : 3,
"K" : 2,
"L" : 6,
"M" : 1,
"N" : 2,
"P" : 4,
"Q" : 2,
"R" : 6,
"S" : 6,
"T" : 4,
"V" : 4,
"W" : 1,
"Y" : 2}

# initialize the number of possible mRNAs with 3 possible STOP codons
mrnas = 3

for aa in protein:
    mrnas *= code_reversed[aa]

print(mrnas)
