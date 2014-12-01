#!/usr/bin/env python3

"""
Expected Number of Restriction Sites
====================================

Given: A positive integer n (n≤1,000,000), a DNA string s of even length
at most 10, and an array A of length at most 20, containing numbers
between 0 and 1.

Return: An array B having the same length as A in which B[i] represents
the expected number of times that s will appear as a substring of a random
DNA string t of length n, where t is formed with GC-content A[i].

Sample Dataset
--------------

10
AG
0.25 0.5 0.75

Sample Output
-------------

0.422 0.563 0.422

"""

import sys

# a function to compute the probability of dna string,
# if we know the gc content (used as the base of probability
# for computing each nucleotide)
def gc_to_dna_probability(gc, dna):
  prob = 1
  for nt in dna:
    if nt in "CG":
      prob *= gc/2
    else:
      prob *= (1-gc)/2
  return prob

# open the file
with open(sys.argv[1], 'r') as in_file:
  n = int(in_file.readline().strip())
  s = in_file.readline().strip().upper()
  a = [float(x) for x in in_file.readline().split()]

# how many ways are there to locate a substring of our
# given length in the long string of the length n?
locations = n + 1 - len(s)

# now compute the numbers of expected occurences
exp_occ = []
for gc in a:
  exp_occ.append(locations * gc_to_dna_probability(gc, s))

# print the results
print(" ".join("{:.4f}".format(x) for x in exp_occ))
