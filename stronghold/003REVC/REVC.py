#!/usr/bin/env python3

"""
Complementing a Strand of DNA
=============================

In DNA strings, symbols 'A' and 'T' are complements of each other, as are 'C'
and 'G'. The reverse complement of a DNA string s is the string sc formed by
reversing the symbols of s, then taking the complement of each symbol (e.g.,
the reverse complement of "GTCA" is "TGAC").

Given: A DNA string s of length at most 1000 bp.

Return: The reverse complement sc of s.

Sample Dataset
--------------

AAAACCCGGT

Sample Output
-------------

ACCGGGTTTT

"""

import sys

with open(sys.argv[1], 'r') as input_file:   # to automagically close the file
  data = input_file.read()                   # when leaving the nested block

dna = ''.join(data.split())   # the split/join trick to remove all whitespace
reverse_original = dna[::-1]  # the extended slice string reversing trick

# complementing, preserving case...
complement_table = str.maketrans("ATCGatcg","TAGCtagc")
reverse_complement = reverse_original.translate(complement_table)

print(reverse_complement),
