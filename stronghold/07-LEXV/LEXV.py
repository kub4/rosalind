#!/usr/bin/env python3

"""
Ordering Strings of Varying Length Lexicographically
====================================================

Given: A permutation of at most 12 symbols defining an ordered alphabet A and a
positive integer n (n≤4).

Return: All strings of length at most n formed from A, ordered
lexicographically. (Note: As in “Enumerating k-mers Lexicographically”,
alphabet order is based on the order in which the symbols are given.)

Sample Dataset
--------------

D N A
3

Sample Output
-------------

D
DD
DDD
DDN
DDA
DN
DND
DNN
DNA
DA
DAD
DAN
DAA
N
ND
NDD
NDN
NDA
NN
NND
NNN
NNA
NA
NAD
NAN
NAA
A
AD
ADD
ADN
ADA
AN
AND
ANN
ANA
AA
AAD
AAN
AAA

"""

import sys, itertools

with open(sys.argv[1], 'r') as in_file:    # to automagically close the file
  alphabet = in_file.readline().split()
  n = int(in_file.readline().strip())

l = []

# create all the strings
for i in range(1,n+1):
  for p in itertools.product(alphabet, repeat=i):
    l.append(''.join(p))

# sorting lambda magic! 
sorted_list = sorted(l, key=lambda w: [alphabet.index(c) for c in w])

# output the results
for i in sorted_list:
  print(i)
