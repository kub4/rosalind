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

def prefix_generator(dna):
  for l in range(len(dna)):
    yield dna[:l+1]  

# open the file and read the sequence
with open(sys.argv[1], 'r') as in_file:
  lines = in_file.readlines()
dna = [x.strip().upper() for x in lines[1:]]
dna = "".join(dna)

"""
# brute force method
failure_array = []
for i, nt in enumerate(dna):
  f = 0
  if i == 0:
    failure_array.append(f) # by convention, first element is zero
  else:
    for l in range(1,i):    # we try all permissible lengths
      if dna[:l] == dna[i-l+1:i+1]:
        if l > f:
          f = l             # and record the longest match
    failure_array.append(f) # then append it to the failure array
"""
"""
# more clever brute force, still hopeless
failure_array = [0]*len(dna)
for prefix in prefix_generator(dna):
  l = len(prefix)»
  for i in range(l,len(dna)):
    if prefix == dna[i-l+1:i+1]:
      failure_array[i] = l
"""

failure_array = [0]*len(dna)

for j in range(1,len (dna)): # first nucleotide to compare with prefix
  if dna[j] == dna [0]:
    failure_array[j] = 1

for j in range(2,len(dna)):
  if failure_array[j-1] == 1 and dna[j] == dna[1]:
    failure_array[j] = 2

for j in range(3,len(dna)):
  if failure_array[j-1] == 2 and dna[j] == dna[2]:
    failure_array[j] = 3

for j in range(4,len(dna)):
  if failure_array[j-1] == 3 and dna[j] == dna[3]:
    failure_array[j] = 4

for j in range(5,len(dna)):
  if failure_array[j-1] == 4 and dna[j] == dna[4]:
    failure_array[j] = 5

for j in range(6,len(dna)):
  if failure_array[j-1] == 5 and dna[j] == dna[5]:
    failure_array[j] = 6


"""
for i in range(1, len(dna)-1): # prefix length
  for j in range(i,len (dna)): # first nucleotide to compare with prefix
    if dna[j:j+i] == dna [0:i]:
      failure_array[j] = i
"""

print(" ".join(map(str,failure_array)))
