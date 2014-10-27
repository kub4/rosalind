#!/usr/bin/env python3

"""
Enumerating Oriented Gene Orderings
===================================

A signed permutation of length n is some ordering of the positive integers
{1,2,…,n} in which each integer is then provided with either a positive or
negative sign (for the sake of simplicity, we omit the positive sign). For
example, π=(5,−3,−2,1,4) is a signed permutation of length 5.

Given: A positive integer n≤6.

Return: The total number of signed permutations of length n, followed by
a list of all such permutations (you may list the signed permutations in
any order).

Sample Dataset
--------------

2

Sample Output
-------------

8
-1 -2
-1 2
1 -2
1 2
-2 -1
-2 1
2 -1
2 1

"""

import sys, math, itertools

with open(sys.argv[1], 'r') as in_file:   # to automagically close the file
  n = int(in_file.read().strip())         # when leaving the nested block

# we need tu multiple ordinary permutations by possible sign arrangements
number_of_permutations = math.factorial(n) * (2 ** n) 

print(number_of_permutations)

for p in itertools.permutations(range(1,n+1)):
  for signs in itertools.product((-1,1),repeat=n):
    print(' '.join(map(str,(p*sign for p, sign in zip(p,signs)))))
