#!/usr/bin/env python3

"""
Perfect Matchings and RNA Secondary Structures
==============================================

Given an RNA string s=s1…sn, a bonding graph for s is formed as follows. First,
assign each symbol of s to a node, and arrange these nodes in order around a
circle, connecting them with edges called adjacency edges. Second, form all
possible edges {A, U} and {C, G}, called basepair edges.

Note that a matching contained in the basepair edges will represent one
possibility for base pairing interactions in s, as shown in Figure 5. For such
a matching to exist, s must have the same number of occurrences of 'A' as 'U'
and the same number of occurrences of 'C' as 'G'.

Given: An RNA string s of length at most 80 bp having the same number of
occurrences of 'A' as 'U' and the same number of occurrences of 'C' as 'G'.

Return: The total possible number of perfect matchings of basepair edges in the
bonding graph of s.

Sample Dataset
--------------

>Rosalind_23
AGCUAGUCAU

Sample Output
-------------

12

"""

import sys
from math import factorial as fac

with open(sys.argv[1], 'r') as in_file:   # to automagically close the file
  data = in_file.read().split()

rna = ''.join(data[1:]).upper()

length = len(rna)
a = rna.count("A")
u = rna.count("U")
c = rna.count("C")
g = rna.count("G")

assert length == a+u+c+g, "Sanity check failed!!"
assert a==u, "Number of A and U should be equal!"
assert c==g, "Number of C and G should be equal!"

"""
and after a little thought...
number of possible perfect mathings is...
"""
print(fac(a)*fac(g))
