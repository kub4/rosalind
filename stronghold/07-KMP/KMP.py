#!/usr/bin/env python3

"""
Speeding Up Motif Finding
=========================

A prefix of a length n string s is a substring s[1:j]; a suffix of s is a
substring s[k:n] (1-based numbering).

The failure array of s is an array P of length n for which P[k] is the length
of the longest substring s[j:k] that is equal to some prefix s[1:k−j+1], where
j cannot equal 1 (otherwise, P[k] would always equal k). By convention,
P[1]=0.

Given: A DNA string s (of length at most 100 kbp) in FASTA format.

Return: The failure array of s.

Sample Dataset
--------------

>Rosalind_87
CAGCATGGTATCACAGCAGAG

Sample Output
-------------

0 0 0 1 2 0 0 0 0 0 0 1 2 1 2 3 4 5 3 0 0

"""

import sys

# open the file and parse the sequence
with open(sys.argv[1], 'r') as in_file:
  lines = in_file.readlines()
dna = [x.strip().upper() for x in lines[1:]]
dna = "".join(dna)

# measure the dna
length = len(dna)

# create the working array
working_array = [[] for i in range(length)]

# here comes the magick
length = len(dna)
for nt in range(1,length):
  if dna[nt] == dna[0]:
    working_array[nt].append(1)
  for x in working_array[nt-1]:
    if dna[nt] == dna[x]:
      working_array[nt].append(x+1)

# create the failure array
failure_array = []
for f in working_array:
  if f:
    failure_array.append(f[-1])
  else:
    failure_array.append(0)

# print the results
print(" ".join(map(str,failure_array)))
