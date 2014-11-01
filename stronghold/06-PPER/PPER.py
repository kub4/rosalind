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
import math

# open the file and extract the numbers
with open(sys.argv[1], 'r') as in_file:
  data = in_file.read().split()

n = int(data[0])
k = int(data[1])

assert (n > 0) and (k > 0) and (n >= k), "Bad input!"

# now to compute the partial permutations, clearly, solution is:
"""
combination (selecting k items from n items)
              *
permutation (of the selected k items)

              |
              v

combination = fac(n)/(fac(k)*fac(n-k))
              *
permutation = fac(k)

              |
              v

        fac(n)/fac(n-k)
"""
# well, that is pretty neat formula and well usable in python 
# which has good bignum support... but to implement the anti-overflow
# modulo trick during the computation, we can go further...
"""
              |
              v
(n-k+1) * ... * (n-2) * (n-1) * n
"""
# so...

pper = 1

for i in range (n-k+1, n+1):
  pper *= i
  pper %= 1000000 # the modulo trick used during computation

print(pper)

"""
as can be tested against the following function,
ẗhe modulo trick makes the computing with large inputs
(like: 800000 2000) much faster even in Python...
"""
#print((math.factorial(n)//math.factorial(n-k))%1000000)
