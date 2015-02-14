#!/usr/bin/env python3

"""
Reverse Complement Problem
==========================

Find the reverse complement of a DNA string.

Sample Dataset
--------------

AAAACCCGGT

Sample Output
-------------

ACCGGGTTTT

"""

import sys

# open the input file and parse the data
with open(sys.argv[1], 'r') as input_file:
  data = input_file.read().split()
  dna = "".join(data)

# and do all the magic in just one line
print(dna[::-1].translate(str.maketrans("ACGT", "TGCA")))
