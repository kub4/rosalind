#!/usr/bin/env python3

"""
Finding a Motif in DNA
======================

Given: Two DNA strings s and t (each of length at most 1 kbp).

Return: All locations of t as a substring of s.
USE 1-BASED NUMBERING!!!!!!!!!!!!!!!!!!!!!!!!!!

Sample Dataset
--------------

GATATATGCATATACTT
ATAT

Sample Output
-------------

2 4 10

"""
import sys

# define a function that finds all occurences of a substring
def find_all_substrings(string, substring):
  start = 0
  while True:
    start = dna.find(substring, start)
    if start == -1: return
    yield start
    start += 1 # we need overlaps

with open(sys.argv[1], 'r') as in_file:   # to automagically close the file
  dna = in_file.readline().strip()        # when leaving the nested block
  motif = in_file.readline().strip()

# create a list of _1-based_ locations
locations = [loc+1 for loc in find_all_substrings(dna, motif)]

print(' '.join(map(str,locations)))
