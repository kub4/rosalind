#!/usr/bin/env python3

"""
Introduction to Set Operations
==============================

Given: A positive integer n (n≤20,000) and two subsets A and B of {1,2,…,n}.

Return: Six sets: A∪B, A∩B, A−B, B−A, Ac, and Bc (where set complements are
taken with respect to {1,2,…,n}).


Sample Dataset
--------------

10
{1, 2, 3, 4, 5}
{2, 8, 5, 10}

Sample Output
-------------

{1, 2, 3, 4, 5, 8, 10}
{2, 5}
{1, 3, 4}
{8, 10}
{8, 9, 10, 6, 7}
{1, 3, 4, 6, 7, 9}

"""

import sys, ast

# open a file and parse the data
# sets are formatted as python set strings,
# so we can use literal_eval to parse them
with open(sys.argv[1], 'r') as in_file:
  data = in_file.read().splitlines()

n = int(data[0])
setA = ast.literal_eval(data[1])
setB = ast.literal_eval(data[2])

setC = set(range(1,n+1))

print(setA | setB)
print(setA & setB)
print(setA - setB)
print(setB - setA)
print(setC - setA)
print(setC - setB)
