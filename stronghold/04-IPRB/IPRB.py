#!/usr/bin/env python3

"""
Mendel's First Law
==================

Given: Three positive integers k, m, and n, representing a population
containing k+m+n organisms: k individuals are homozygous dominant for a factor,
m are heterozygous, and n are homozygous recessive.

Return: The probability that two randomly selected mating organisms will
produce an individual possessing a dominant allele (and thus displaying the
dominant phenotype). Assume that any two organisms can mate.

Sample Dataset
--------------

2 2 2

Sample Output
-------------

0.78333

"""
import sys

with open(sys.argv[1], 'r') as in_file:   # to automagically close the file
  data = list(map(int,in_file.readline().split()))

k = data[0]   # dominant homozygotes
m = data[1]   # heterozygotes
n = data[2]   # recessive homozygotes
s = k+m+n     # sum of all individuals

assert (k > 0), "k should be positive!"
assert (m > 0), "m should be positive!"
assert (n > 0), "n should be positive!"
assert (s > 1), "need a pair to have sex!"

dom_child_probability = 0 # we start with zero

# if "father" belongs to k, "mother" does not matter
dom_child_probability += (k/s * 1)

# if "father" belongs to n, "mother" decides
dom_child_probability += (n/s * ((k/(s-1)) + (0.5*m)/(s-1)))

# if "father" belongs to m, both parents are important
### (1) "father" belongs to m, "mother" to k:
dom_child_probability += (m/s * k/(s-1))
### (2) "father" belongs to m, "mother" to n:
dom_child_probability += 0.5 * (m/s * n/(s-1))
### (2) "father" belongs to m, "mother" also to m:
dom_child_probability += 0.75 * (m/s * (m-1)/(s-1))

print(dom_child_probability)
