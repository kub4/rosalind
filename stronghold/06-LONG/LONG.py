#!/usr/bin/env python3

"""
Genome Assembly as Shortest Superstring
=======================================

For a collection of strings, a larger string containing every one of the
smaller strings as a substring is called a superstring. By the assumption of
parsimony, a shortest possible superstring over a collection of reads serves
as a candidate chromosome.

Given: At most 50 DNA strings whose length does not exceed 1 kbp in FASTA
format (which represent reads deriving from the same strand of a single linear
chromosome).

The dataset is guaranteed to satisfy the following condition: there exists a
unique way to reconstruct the entire chromosome from these reads by gluing
together pairs of reads that overlap by more than half their length.

Return: A shortest superstring containing all the given strings (thus
corresponding to a reconstructed chromosome).

Sample Dataset
--------------

>Rosalind_56
ATTAGACCTG
>Rosalind_57
CCTGCCGGAA
>Rosalind_58
AGACCTGCCG
>Rosalind_59
GCCGGAATAC

Sample Output
-------------

ATTAGACCTGCCGGAATAC

"""

import sys

# open a file and parse the data
with open(sys.argv[1], 'r') as in_file:
  lines = in_file.read().splitlines()

sequences=[]
for line in lines:
  if len(line) > 0 and line[0]==">":
    sequences.append([])
  elif sequences:
    sequences[-1].append(line)
sequences = ["".join(s).upper() for s in sequences]

# create jointable
jointable = []
for i, s in enumerate(sequences):
  l = (len(s)+1)//2
  for j, t in enumerate(sequences):
    if i != j:
      match = t.find(s[-l:])
      if match >= 0 and s[-l-match:]==t[:l+match]:
        jointable.append((i,j,l+match)) 

left_fragments = [x[0] for x in jointable]
right_fragments = [x[1] for x in jointable]
start = set(left_fragments) - set(right_fragments)
assert len(left_fragments) == len(right_fragments), "This is guaranteed not to happen!"
assert len(start)==1, "This is guaranteed not to happen!"

chromosome = ""
fragment_to_add = list(start)[0]
while True:
  chromosome += sequences[fragment_to_add]
  x = [x for x in jointable if x[0]==fragment_to_add]
  if x:
    fragment_to_add = x[0][1]
    chromosome = chromosome[:-x[0][2]]
  else:
    break

print(chromosome)
