#!/usr/bin/env python3

"""
Majority Element
================

An array A[1..n] is said to have a majority element if more
than half of its entries are the same.

Given: A positive integer k≤20, a positive integer n≤104,
and k arrays of size n containing positive integers not
exceeding 105.

Return: For each array, output an element of this array
occurring strictly more than n/2 times if such element exists,
and "-1" otherwise.

Sample Dataset
--------------

4 8
5 5 5 5 5 5 5 5
8 7 7 7 1 7 3 7
7 1 6 5 10 100 1000 1
5 1 6 7 1 1 10 1

Sample Output
-------------

5 7 -1 -1

"""
import sys
from collections import Counter

# open the input file and parse the data
with open(sys.argv[1], 'r') as in_file:
  k,n = [int(x) for x in in_file.readline().split()]
  arrays = []
  for a in range(k):
    arrays.append([int(x) for x in in_file.readline().split()])

# minimum number of identical elements to be called majority
min_number = n//2 + 1

# create the majority elemens array, initialized with -1's
majority_elements = [-1]*k

# with python, it is easy to cheat :)
for i,a in enumerate(arrays):
  most_common = Counter(a).most_common(1)[0]
  if most_common[1] >= min_number:
    majority_elements[i] = most_common[0]

# print out the results
print(" ".join(map(str,majority_elements)))
