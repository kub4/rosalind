#!/usr/bin/env python3

"""
Partial Permutations
====================

A partial permutation is an ordering of only k objects taken from a collection
containing n objects (k≤n). For example, one partial permutation of three
of the first eight positive integers is given by (5,7,2).

The statistic P(n,k) counts the total number of partial permutations of
k objects that can be formed from a collection of n objects. Note that
P(n,n) is just the number of permutations of n objects, which is equal to n!.

Given: Positive integers n and k such that 100≥n>0 and 10≥k>0.

Return: The total number of partial permutations P(n,k), modulo 1,000,000.

Sample Dataset
--------------

21 7

Sample Output
-------------

51200

"""

import sys
from math import factorial as fac

# open the file and extract the numbers
with open(sys.argv[1], 'r') as in_file:
  data = in_file.read().split()

n = int(data[0])
k = int(data[1])

assert (n > 0) and (k > 0) and (n >= k), "Bad input!"

# now to compute the partial permutations
# clearly, solution is:
"""
combination (selecting k items from n items)
*
permutation (of the selected k items)
"""
combination = int(fac(n)/(fac(k)*fac(n-k)))
permutation = fac(k)

print(combination * permutation % 1000000)
