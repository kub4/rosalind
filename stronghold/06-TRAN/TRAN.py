#!/usr/bin/env python3

"""
Transitions and Transversions
=============================

For DNA strings s1 and s2 having the same length, their transition/transversion
ratio R(s1,s2) is the ratio of the total number of transitions to the total
number of transversions, where symbol substitutions are inferred from
mismatched corresponding symbols as when calculating Hamming distance.

Given: Two DNA strings s1 and s2 of equal length (at most 1 kbp).

Return: The transition/transversion ratio R(s1,s2).

Sample Dataset
--------------

>Rosalind_0209
GCAACGCACAACGAAAACCCTTAGGGACTGGATTATTTCGTGATCGTTGTAGTTATTGGA
AGTACGGGCATCAACCCAGTT
>Rosalind_2200
TTATCTGACAAAGAAAGCCGTCAACGGCTGGATAATTTCGCGATCGTGCTGGTTACTGGC
GGTACGAGTGTTCCTTTGGGT

Sample Output
-------------

1.21428571429

"""

import sys

# open the file and create a list of sequences
with open(sys.argv[1], 'r') as in_file:
  data = in_file.read().split(">")
  dnas = ["".join(x.upper().split()[1:]) for x in filter(None, data)]

# measure length, sanity check
l = len(dnas[0])
assert len(dnas[1]) == l, "Both strings must have the same length!"

# count transitions and transversions
transitions = 0
transversions = 0
purines = ["A","G"]

for i in range(l):
  if dnas[0][i] == dnas[1][i]:
    continue
  elif (dnas[0][i] in purines) == (dnas[1][i] in purines):
    transitions += 1
  else:
    transversions += 1

# do not divide by zero
if transversions > 0:
  r = transitions/transversions
else:
  r = "{0} transitions and no transversions!".format(transitions)

# print the results
print(r)
