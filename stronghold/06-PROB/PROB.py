#!/usr/bin/env python3

"""
Introduction to Random Strings
==============================

Given: A DNA string s of length at most 100 bp and an array A containing at
most 20 numbers between 0 and 1 representing various GC-contents.

Return: An array B having the same length as A in which B[k] represents the
common logarithm of the probability that a random string constructed with the
GC-content found in A[k] will match s exactly.

Sample Dataset
--------------

ACGATACAA
0.129 0.287 0.423 0.476 0.641 0.742 0.783

Sample Output
-------------

-5.737 -5.217 -5.263 -5.360 -5.958 -6.628 -7.009

"""

import sys
import math

# open the file and extract the data
with open(sys.argv[1], 'r') as in_file:
  dna = in_file.readline().upper().strip()
  gcc = [float(x) for x in in_file.readline().split()]

# make some sanity checks
for gc in gcc:
  assert 0 < gc < 1, "Need 1.0 > GC > 0.0 !"
  # for gc==1 or gc==0, special handling would be needed because of the use
  # of the logarithm function and we would need to output something like minus
  # infinity for zero probabilities, only GC-only or AT-only strings with
  # matching GC content of 1.0 or 0.0 would have probability equal to 1...
  # ...but it does not make much biological sense and rosalind does not seem
  # to use gc==1 or gc==0 in the inputs anyway...
assert set(dna) <= set("ACGT"), "Not a valid DNA string!"

"""
# compute the probabilities normally, log them lately - this is bad for long
# (==improbable) strings due to limited resolution of floats
# uncomment this for comparision
probs_lin = []
for gc in gcc:
  prob = 1 # initialize string probability (for multiplication)
  for nt in dna:
    if nt in ["G","C"]:
      prob *= gc/2
    else:
      prob *= (1-gc)/2
  probs_lin.append(prob)

print(" ".join(map(str,probs_lin)))
print(" ".join(map(str,[math.log(x,10) if x > 0 else None for x in probs_lin])))
# with a random string of 500 nucleotides, I got zero/None for all gc==0.25
# and lower...
"""

# compute the probabilities as a sum of logarithms, usable
# even for long strings
probs_log = []
for gc in gcc:
  prob = 0 # initialize string probability (for addition)
  for nt in dna:
    if nt in ["G","C"]:
      prob += math.log(gc/2,10)
    else:
      prob += math.log((1-gc)/2,10)
  probs_log.append(prob)

print(" ".join(map(str,probs_log)))
