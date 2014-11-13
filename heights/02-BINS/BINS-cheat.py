#!/usr/bin/env python3

"""
Binary Search
=============

Given: Two positive integers n≤10^5 and m≤10^5, a sorted array A[1..n] of
integers from −10^5 to 10^5 and a list of m integers −10^5≤k1,k2,…,km≤10^5.

Return: For each ki, output an index 1≤j≤n s.t. A[j]=ki
or "-1" if there is no such index.

!!!1-based numbering!!!

Sample Dataset
--------------

5
6
10 20 30 40 50
40 10 35 15 40 20

Sample Output
-------------

4 1 -1 -1 4 2

"""
import sys

# open the input file and parse the data
with open(sys.argv[1], 'r') as in_file:
  n = int(in_file.readline().strip())
  m = int(in_file.readline().strip())
  array = [int(x) for x in in_file.readline().split()]
  klist = [int(x) for x in in_file.readline().split()]

# sanity check
assert n == len(array) and m == len(klist), "Bad input!" 

# with python, it is easy to cheat :)
indexes = []
for k in klist:
  if k in array:
    indexes.append((array.index(k)+1)) # 1-based numbering
  else:
    indexes.append(-1)

# print out the results
print(" ".join(map(str,indexes)))
