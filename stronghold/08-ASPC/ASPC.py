#!/usr/bin/env python3

"""
Introduction to Alternative Splicing
====================================

Given: Positive integers n and m with 0≤m≤n≤2000.

Return: The sum of combinations C(n,k) for all k satisfying m≤k≤n,
modulo 1,000,000.

Sample Dataset
--------------

6 3

Sample Output
-------------

42

"""

import sys

# open the file
with open(sys.argv[1], 'r') as in_file:
  n, m = [int(x) for x in in_file.readline().split()]

"""
so the unoptimized formula for combinations is:
            n!
C(n,k) = ---------
         k!*(n-k)!
"""
print(n,m)
