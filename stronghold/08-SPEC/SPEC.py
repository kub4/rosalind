#!/usr/bin/env python3

"""
Inferring Protein from Spectrum
===============================

The prefix spectrum of a weighted string is the collection
of all its prefix weights.

Given: A list L of n (n≤100) positive real numbers.

Return: A protein string of length n−1 whose prefix spectrum
is equal to L (if multiple solutions exist, you may output
any one of them). Consult the monoisotopic mass table.

Sample Dataset
--------------

3524.8542
3710.9335
3841.974
3970.0326
4057.0646

Sample Output
-------------

WMQS

"""

import sys, random

# the monoisotopic mass table
monoisotopic_masses = {
  "A" :  71.03711,
  "C" : 103.00919,
  "D" : 115.02694,
  "E" : 129.04259,
  "F" : 147.06841,
  "G" :  57.02146,
  "H" : 137.05891,
  "I" : 113.08406,
  "K" : 128.09496,
  "L" : 113.08406,
  "M" : 131.04049,
  "N" : 114.04293,
  "P" :  97.05276,
  "Q" : 128.05858,
  "R" : 156.10111,
  "S" :  87.03203,
  "T" : 101.04768,
  "V" :  99.06841,
  "W" : 186.07931,
  "Y" : 163.06333    }

# get the spectrum from the input file
spectrum = []
with open(sys.argv[1], 'r') as in_file:
  for line in in_file:
    if line.strip():
      spectrum.append(float(line.strip()))

# sort the spectrum
spectrum = sorted(spectrum)

# create a list for the protein
protein = []

# analyze the spectrum
for i in range(1, len(spectrum)):
  w = spectrum[i]-spectrum[i-1]
  valid_aas = []
  # compare rounded numbers to get the matches
  for aa, mass in monoisotopic_masses.items():
    if round(mass,2) == round(w,2):
      valid_aas.append(aa)
  #for more fun, if more choices, choose randomly
  protein.append(random.choice(valid_aas))

# print the resulting protein
print("".join(protein))
