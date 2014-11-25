#!/usr/bin/env python3

"""
Matching Random Motifs
======================

Given: A positive integer N≤100000, a number x between 0 and 1, and a DNA
string s of length at most 10 bp.

Return: The probability that if N random DNA strings having the same length as
s are constructed with GC-content x, then at least one of the strings equals s.
We allow for the same random string to be created more than once. A random
string is constructed so that the probability of choosing each subsequent
symbol is based on a fixed underlying symbol frequency.

Sample Dataset
--------------

90000 0.6
ATAGCCGA

Sample Output
-------------

0.689

"""

import sys

# read data from the input file
with open(sys.argv[1], 'r') as in_file:
  numbers = in_file.readline().split()
  dna = in_file.readline().strip().upper()

# parse numbers
N = int(numbers[0])
x = float(numbers[1])

# compute the probability of one string matching s
prob_one = 1
for nt in dna:
  if nt in ["G", "C"]:
    prob_one *= x/2
  elif nt in ["A", "T"]:
    prob_one *= (1-x)/2
  else:
    raise Exception("Bad DNA!")

# compute the probability of one string not matching s
neg_prob_one = 1 - prob_one

# compute the probabily of none of N strings matching s
neg_prob_all = neg_prob_one**N

# compute the probability of at least one of N matching s
prob_at_least_one = 1-neg_prob_all

# print the result with 3 decimal places
print("{:.3f}".format(prob_at_least_one))
