#!/usr/bin/env python3

"""
Counting DNA Nucleotides
========================

A string is simply an ordered collection of symbols selected from some alphabet
and formed into a word; the length of a string is the number of symbols that it
contains. An example of a length 21 DNA string (whose alphabet contains the
symbols 'A', 'C', 'G', and 'T') is "ATGCTTCAGAAAGGTCTTACG."

Given: A DNA string s of length at most 1000 nt.

Return: Four integers (separated by spaces) counting the respective number of
times that the symbols 'A', 'C', 'G', and 'T' occur in s.

"""

import sys

with open(sys.argv[1], 'r') as input_file:   # to automagically close the file
  data = input_file.read()                   # when leaving the nested block

data = ''.join(data.split()) # the split/join trick to remove all whitespace
data = data.upper()          # convert to uppercase to simplify counting  

for c in ['A','C','G']:
    print(data.count(c), end=" ")   # separate by spaces, no newlines
print(data.count('T'))              # finish with a newline
