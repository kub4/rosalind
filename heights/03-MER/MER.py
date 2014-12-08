#!/usr/bin/env python3

"""
================

The merging procedure is an essential part of “Merge Sort”
(which is considered in one of the next problems).

Given: A positive integer n≤10^5 and a sorted array A[1..n] of integers
from −10^5 to 10^5, a positive integer m≤10^5 and a sorted array B[1..m]
of integers from −10^5 to 10^5.

Return: A sorted array C[1..n+m] containing all the elements of A and B.

Sample Dataset
--------------

4
2 4 10 18
3
-5 11 12

Sample Output
-------------

-5 2 4 10 11 12 18

"""
import sys

# open the input file and parse the data
with open(sys.argv[1], 'r') as in_file:
  lines = in_file.readlines()

a = [int(x) for x in lines[1].split()]
b = [int(x) for x in lines[3].split()]

# python is python .)
print(" ".join(str(x) for x in sorted(a+b)))
