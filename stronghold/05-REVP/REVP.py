#!/usr/bin/env python3

"""
Locating Restriction Sites
==========================

A DNA string is a reverse palindrome if it is equal to its reverse complement.
For instance, GCATGC is a reverse palindrome because its reverse complement is
GCATGC.

Given: A DNA string of length at most 1 kbp in FASTA format.

Return: The position and length of every reverse palindrome in the string
having length between 4 and 12. You may return these pairs in any order.

Sample Dataset
--------------

>Rosalind_24
TCAATGCATGCGGGTCTATATGCAT

Sample Output
-------------

4 6
5 4
6 6
7 4
17 4
18 4
20 6
21 4

"""
import sys

# extract dna from the file
with open(sys.argv[1], 'r') as in_file:
  dna = ''.join(in_file.read().upper().split("\n")[1:])

# create the reverse complement dna string
basepairs = str.maketrans("ACGT", "TGCA")
rev_dna = dna[::-1].translate(basepairs)

# measure the dna strand
dl = len(dna)

# search for the reverse palindromes
for pl in range(4,13,2): # reverse palindrom length must be even!
  for i in range(dl-pl+1): # end just in time before the strand ends
    if dna[i:i+pl] == rev_dna[dl-i-pl:dl-i]:
      print(i+1, pl) # print 1-based numbering
