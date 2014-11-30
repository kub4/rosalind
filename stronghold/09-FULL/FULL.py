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

# the inverted monoisotopic mass table, because I and L
# have the same mass, I should be treated as I or L if
# it is important
monoisotopic_masses = {
  57.02146 : "G",
  71.03711 : "A",
  87.03203 : "S",
  97.05276 : "P",
  99.06841 : "V",
 101.04768 : "T",
 103.00919 : "C",
 113.08406 : "I", # or L
 114.04293 : "N",
 115.02694 : "D",
 128.05858 : "Q",
 128.09496 : "K",
 129.04259 : "E",
 131.04049 : "M",
 137.05891 : "H",
 147.06841 : "F",
 156.10111 : "R",
 163.06333 : "Y",
 186.07931 : "W"      }

def almostequal(number1, number2):
  delta = 0.01
  if abs(number1 - number2) < delta:
    return True
  else:
    return False

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
total_mass = masses[0]

# sort the remaining mass list
fragment_masses = sorted(masses[1:])

# select the shortest mass as the bare prefix (may be it is a suffix,
# but we cannot distinguish it
bare_prefix = fragment_masses.pop(0)

#print(fragment_masses)
for ii,x in enumerate(fragment_masses):
  #print(",?")
  if almostequal(total_mass-bare_prefix, x):
    #print(total_mass-bare_prefix, x)
    fragment_masses[ii] = False
    #rint("pl\n\npppppp!!!!!!!!!!")
    break

protein = []
curr_mass = bare_prefix
#print(fragment_masses)

while True:
  for i, m in enumerate(fragment_masses):
    #print("outerforbegins")
    #print(fragment_masses)
    #print(curr_mass)
    for k, v in monoisotopic_masses.items():
      if almostequal((m - curr_mass),k):
        #print(m-curr_mass, k, "equal?")
        protein.append(v)
        #print(protein)
        curr_mass = m
        fragment_masses[i] = False
        for y,x in enumerate(fragment_masses):
          if almostequal(total_mass-m, x):
            fragment_masses[y] = False
            break
        print(fragment_masses)
        break
  break
#print(fragment_masses)
assert n == len(protein)
print("".join(protein))

