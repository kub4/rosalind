#!/usr/bin/env python3

"""
Enumerating Gene Orders
=======================

A permutation of length n is an ordering of the positive integers {1,2,…,n}.
For example, π=(5,3,2,1,4) is a permutation of length 5.

Given: A positive integer n≤7.

Return: The total number of permutations of length n, followed by a list
of all such permutations (in any order).

Sample Dataset
--------------

3

Sample Output
-------------

6
1 2 3
1 3 2
2 1 3
2 3 1
3 1 2
3 2 1

"""
import sys, math, itertools

with open(sys.argv[1], 'r') as in_file:   # to automagically close the file
  n = int(in_file.read().strip())         # when leaving the nested block

number_of_permutations = math.factorial(n)

print(number_of_permutations)

for p in itertools.permutations(range(1,n+1)):
  i = 0
  for r in p:
    i += 1
    if i < n: print(r, end=" ")
    else: print(r)
