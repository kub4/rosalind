#!/usr/bin/env python3

"""
Counting Subsets
================

Given: A positive integer n (n≤1000).

Return: The total number of subsets of {1,2,…,n} modulo 1,000,000.
(Including the empty set and the set itself.)


Sample Dataset
--------------

3

Sample Output
-------------

8

"""

import sys
#from math import factorial as fac

# open a file and get n
with open(sys.argv[1], 'r') as in_file:
  n = int(in_file.read().strip())

# this works, but...
"""
# count the subsets
subsets = 0
for k in range(0, n+1):
  subsets += (fac(n)//(fac(k)*fac(n-k)))

# in Python, we can safely modulo in the very end
print(subsets % 1000000)
"""

# the mathematical truth is that the number
# of subsets is equal to n!... and the pow
# function in python supports modulo...
print(pow(2,n,1000000))
