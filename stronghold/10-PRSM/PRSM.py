#!/usr/bin/env python3

"""
Matching a Spectrum to a Protein
================================

The complete spectrum of a weighted string s is the multiset S[s]
containing the weights of every prefix and suffix of s.

Given: A positive integer n followed by a collection of n protein strings s1,
s2, ..., sn and a multiset R of positive numbers (corresponding to the complete
spectrum of some unknown protein string).

Return: The maximum multiplicity of R⊖S[sk] taken over all strings sk, followed
by the string sk for which this maximum multiplicity occurs (you may output any
such value if multiple solutions exist).

Sample Dataset
--------------

4
GSDMQS
VWICN
IASWMQS
PVSMGAD
445.17838
115.02694
186.07931
314.13789
317.1198
215.09061

Sample Output
-------------


"""

import sys
from collections import Counter

# the monoisotopic mass table
monoisotopic_masses = {
"A" : 71.03711,
"C" : 103.00919,
"D" : 115.02694,
"E" : 129.04259,
"F" : 147.06841,
"G" : 57.02146,
"H" : 137.05891,
"I" : 113.08406,
"K" : 128.09496,
"L" : 113.08406,
"M" : 131.04049,
"N" : 114.04293,
"P" : 97.05276,
"Q" : 128.05858,
"R" : 156.10111,
"S" : 87.03203,
"T" : 101.04768,
"V" : 99.06841,
"W" : 186.07931,
"Y" : 163.06333 }


rnd = 3
def weight(p):
  w = 0
  for aminoacid in p:
    w += int(round(100000*monoisotopic_masses[aminoacid])+0.1)
  return w

def complete_spectrum(p):
  l = len(p)
  spectrum = []
  for prefix in range(1,l+1): # l if complete string is not its own substring :)
    spectrum.append(weight(p[:prefix]))
  for suffix in range(0,l): # 1 if complete string is not its own substring :)
    spectrum.append(weight(p[suffix:]))
  return spectrum
  

# open and read the file
with open(sys.argv[1], 'r') as in_file:
  lines = in_file.readlines()

# parse the data
n        = int(lines[0].strip())
peptides = [p.strip() for p in lines[1:n+1]]
weights  = [int(round(float(w)*100000)+0.1) for w in lines[n+1:]]


max_mult = []
for p in peptides:
  difference = []
  for w in weights:
    for v in complete_spectrum(p):
      difference.append(w-v)
  c=Counter(difference)
  max_mult.append(c.most_common(1)[0][1])

print(max(max_mult))
print(peptides[max_mult.index(max(max_mult))])
