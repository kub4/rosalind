#!/usr/bin/env python3

"""
Enumerating k-mers Lexicographically
====================================

Given: A collection of at most 10 symbols defining an ordered alphabet, and a
positive integer n (n≤10).

Return: All strings of length n that can be formed from the alphabet, ordered
lexicographically.

Sample Dataset
--------------

T A G C
2

Sample Output
-------------

TT
TA
TG
TC
AT
AA
AG
AC
GT
GA
GG
GC
CT
CA
CG
CC

"""

import sys, itertools

with open(sys.argv[1], 'r') as in_file:    # to automagically close the file
  alphabet = in_file.readline().split()
  n = int(in_file.readline().strip())

for p in itertools.product(alphabet, repeat=n):
  print("".join(p))
