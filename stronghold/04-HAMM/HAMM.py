#!/usr/bin/env python3

"""
Counting Point Mutations
========================

Given two strings s and t of equal length, the Hamming distance between s and
t, denoted dH(s,t), is the number of corresponding symbols that differ in s
and t.

Given: Two DNA strings s and t of equal length (not exceeding 1 kbp).

Return: The Hamming distance dH(s,t).

Sample Dataset
--------------

GAGCCTACTAACGGGAT
CATCGTAATGACGGCCT

Sample Output
-------------

7

"""

import sys

with open(sys.argv[1], 'r') as input_file:   # to automagically close the file
  dna1 = input_file.readline().strip()       # when leaving the nested block
  dna2 = input_file.readline().strip()

assert len(dna1) == len(dna2)

hammdist = 0

for base1, base2 in zip(dna1, dna2):
  if (base1 != base2) : hammdist += 1

print(hammdist)
