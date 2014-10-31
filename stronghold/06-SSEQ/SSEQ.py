#!/usr/bin/env python3

"""
Finding a Spliced Motif
=======================

A subsequence of a string is a collection of symbols contained in order
(though not necessarily contiguously) in the string (e.g., ACG is a subsequence
of TATGCTAAGATC). The indices of a subsequence are the positions in the string
at which the symbols of the subsequence appear; thus, the indices of ACG in
TATGCTAAGATC can be represented by (2, 5, 9) - USE 1-BASED COUNTING!!!

As a substring can have multiple locations, a subsequence can have multiple
collections of indices, and the same index can be reused in more than one
appearance of the subsequence; for example, ACG is a subsequence of AACCGGTT
in 8 different ways.

Given: Two DNA strings s and t (each of length at most 1 kbp) in FASTA format.

Return: One collection of indices of s in which the symbols of t appear
as a subsequence of s. If multiple solutions exist, you may return any one.

Sample Dataset
--------------

>Rosalind_14
ACGTACGTGACG
>Rosalind_18
GTA

Sample Output
-------------

3 8 10

"""

import sys

# open the file and create a list of sequences
with open(sys.argv[1], 'r') as in_file:
  data = in_file.read().split(">")
  dnas = ["".join(x.upper().split()[1:]) for x in filter(None, data)]

# measure strands, sanity check
len_dna = len(dnas[0])
len_sub = len(dnas[1])
assert len_dna > len_sub, "The first string should be longer!"

# find the subsequence indices
indices = []

start = 0
for nt in dnas[1]:
  for i in range(start,len_dna):
    if nt == dnas[0][i]:
      indices.append(i+1) # 1-based numbering
      start = i + 1
      break

# check the result, print
if len(indices) == len_sub:
  print(" ".join(map(str, indices)))
else:
  print("Not found")
