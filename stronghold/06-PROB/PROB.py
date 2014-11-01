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
  dna = in_file.readline().strip().upper()
  gcc = [float(x) for x in in_file.readline().split()]

for gc in gcc:
  assert (gc <=1) and (gc >= 0), "Invalid GC content!"


probs_lin = []
for gc in gcc:
  prob = 1 # initialize string probability
  for nt in dna:
    if nt in ["G","C"]:
      prob *= gc/2
    else:
      prob *= (1-gc)/2
  probs_lin.append(prob)

print(probs_lin)
print([math.log(x,10) if x> 0 else None for x in probs_lin])

probs_log = []
for gc in gcc:
  prob = 0 # initialize string probability
  for nt in dna:
    if nt in ["G","C"]:
      prob += math.log(gc/2,10)
    else:
      prob += math.log((1-gc)/2,10)
  probs_log.append(prob)

print(probs_log)
