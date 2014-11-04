#!/usr/bin/env python3
"""
Generates a random string of DNA with equal probability
of A, C, G, or T for each position of a given length.
The output is one line, no FASTA header.
"""

import sys, random

length = int(sys.argv[1])
assert length > 0, "Need length > 0"
dna = []

for nt in range(length):
  dna.append(random.choice("ACGT"))

print("".join(dna))  
