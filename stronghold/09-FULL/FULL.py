#!/usr/bin/env python3

"""
Inferring Peptide from Full Spectrum
====================================

Say that we have a string s containing t as an internal substring, so that
there exist nonempty substrings s1 and s2 of s such that s can be written as
s1ts2. A t-prefix contains all of s1 and none of s2; likewise, a t-suffix
contains all of s2 and none of s1.

Given: A list L containing 2n+3 positive real numbers (n≤100).  The first
number in L is the parent mass of a peptide P, and all other numbers represent
the masses of some b-ions and y-ions of P (in no particular order). You may
assume that if the mass of a b-ion is present, then so is that of its
complementary y-ion, and vice-versa.

Return: A protein string t of length n for which there exist two positive real
numbers w1 and w2 such that for every prefix p and suffix s of t, each of
w(p)+w1 and w(s)+w2 is equal to an element of L. (In other words, there exists
a protein string whose t-prefix and t-suffix weights correspond to the
non-parent mass values of L.) If multiple solutions exist, you may output any
one.

Sample Dataset
--------------

1988.21104821
610.391039105
738.485999105
766.492149105
863.544909105
867.528589105
992.587499105
995.623549105
1120.6824591
1124.6661391
1221.7188991
1249.7250491
1377.8200091

Sample Output
-------------

KEKEP

"""

import sys, random

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

# define the function for analyzing prefix spectra
def prefix_spectrum_to_protein(spectrum):
  print(spectrum)
  """
  From a given, complete, sorted and unpolluted prefix spectrum containing
  zero to all aminoacids of the wanted protein, returns the protein sequence.
  Uses random choice where more than one aminoacid is possible (I/L).
  Requires the monoisotopic_masses dictionary (aminoacids as keys).
  """
  protein = []
  for i in range(1, len(spectrum)):
    w = spectrum[i]-spectrum[i-1]
    valid_aas = []
    print(w)
    # compare rounded numbers to get the matches
    for aa, mass in monoisotopic_masses.items():
      if round(mass,2) == round(w,2):
        valid_aas.append(aa)
    #for more fun, if more choices, choose randomly
    if valid_aas:
      protein.append(random.choice(valid_aas))
    else:
      protein.append("X")
  return(protein)

# create a list to hold all the masses
masses = []

# open the file
with open(sys.argv[1], 'r') as in_file:
  for line in in_file:
    if line.strip():
      masses.append(float(line.strip()))

# compute the length of our wanted subprotein t (which is equal to n)
# (based on the given list description in the problem definition, or,
# in other words, we have all the substrings and no chaos mixed in)
n = (len(masses)-3)/2
assert n==int(n), "input error, bad length of the masses list"
n = int(n)

# separate the mass of the whole s protein from the masses list
smass = masses.pop(0)

# sort the remaining mass list
masses = sorted(masses)

# we have redundant information, so we can split the list in two
prefix_spectrum = masses[:n+1]
suffix_spectrum = list(reversed(masses[n+1:]))

# now build the t subprotein from the both spectra
# but if using a suffix spectrum, it needs to be reversed
t1 = prefix_spectrum_to_protein(prefix_spectrum)
t2 = list(reversed(prefix_spectrum_to_protein(suffix_spectrum)))

# t1 and t2 should be identical, of course
#assert t1==t2, "mismatch of proteins inferred from prefix and suffix spectra"

# print the result
print(t1)
print(t2)
print("".join(t1))
