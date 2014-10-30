#!/usr/bin/env python3

"""
Consensus and Profile
=====================

Say that we have a collection of DNA strings, all having the same length n.
Their profile matrix is a 4×n matrix P in which P1,j represents the number of
times that 'A' occurs in the j-th position of one of the strings, P2,j
represents the number of times that C occurs in the j-th position, and so on.

A consensus string c is a string of length n formed from our collection by
taking the most common symbol at each position; the j-th symbol of c therefore
corresponds to the symbol having the maximum value in the j-th column of the
profile matrix. Of course, there may be more than one most common symbol,
leading to multiple possible consensus strings.

Given: A collection of at most 10 DNA strings of equal length (at most 1 kbp)
in FASTA format.

Return: A consensus string and profile matrix for the collection. (If several
possible consensus strings exist, then you may return any one of them.)

Sample Dataset
--------------

>Rosalind_1
ATCCAGCT
>Rosalind_2
GGGCAACT
>Rosalind_3
ATGGATCT
>Rosalind_4
AAGCAACC
>Rosalind_5
TTGGAACT
>Rosalind_6
ATGCCATT
>Rosalind_7
ATGGCACT


Sample Output
-------------

ATGCAACT
A: 5 1 0 0 5 5 0 0
C: 0 0 1 4 2 0 6 1
G: 1 1 6 3 0 1 0 0
T: 1 5 0 0 0 1 1 6

"""
import sys
from collections import Counter

# open the file and extract useful data
with open(sys.argv[1], 'r') as in_file:
  data = in_file.read().upper().strip()

raw_seq = filter(None,data.split(">"))    # filter removes empty strings
sequences = []                            # a list of lists [[name, sequence]]

for s in raw_seq:
  s = s.split("\n")
  sequences.append([s[0],''.join(s[1:])])

# some sanity checking
length = len(sequences[0][1])
for s in sequences:
  assert len(s[1]) == length, "All DNA strings must have the same lenght!"

# define base order
a = "ACGT"

# create a profile matrix
profile = []
for i in range(length):
  b = [] # bases in the i-th position of all strings
  for s in sequences:
    b.append(s[1][i])
  b = Counter(b)
  counts = (b[a[0]], b[a[1]], b[a[2]], b[a[3]])
  profile.append(counts) # appending count for the i-th position

# get the consensus string
consensus = []
for i in range(length):
  max_val = max(profile[i])
  max_idx = profile[i].index(max_val)
  consensus.append(a[max_idx])
consensus = "".join(consensus)

# printing it all
print(consensus)
for i in range(4):
  print(a[i], end=": ")
  row = []
  for j in range(length):
    row.append(str(profile[j][i]))
  print(" ".join(row))
