#!/usr/bin/env python3

"""
Independent Alleles
===================

Given: Two positive integers k (k≤7) and N (N≤2^k). In this problem, we begin
with Tom, who in the 0th generation has genotype Aa Bb. Tom has two children in
the 1st generation, each of whom has two children, and so on. Each organism
always mates with an organism having genotype Aa Bb.

Return: The probability that at least N Aa Bb organisms will belong to the k-th
generation of Tom's family tree (don't count the Aa Bb mates at each level).
Assume that Mendel's second law holds for the factors.

Sample Dataset
--------------

2 1

Sample Output
-------------

0.684

"""
import sys
from math import factorial as fac

with open(sys.argv[1], 'r') as in_file:   # to automagically close the file
  data = list(map(int,in_file.readline().split()))

k = data[0]
N = data[1]

assert (k>0 and N>0), "k and N must be positive!"

# first, count the children in the k-th generation
n = 2**k

# breeding with Aa Bb results always in 25% Aa Bb
p = 0.25

'''
now we use http://en.wikipedia.org/wiki/Binomial_distribution
cumulative distribution function, for max 'g' Aa Bb organisms,
where g = N-1,.. we will then subtract this probability from 1
'''
#instead of the wiki's 'k', we use g; while n and p are the same

g = N-1

x = 0 # the initialization of probability for max 'g' Aa Bb organisms

for i in range(0,g+1):
  x += (fac(n)/fac(i)/fac(n-i))*(p**i)*((1-p)**(n-i))

print(1-x)

'''
this works for the purposes of the exercise, however
the precission is low, and even within the predetermined
limits it is possible to get for example negative probabilities!
for example for '7 80' I get: -4.440892098500626e-16
'''
