#!/usr/bin/env python3

"""
Calculating Expected Offspring
==============================

Given: Six positive integers, each of which does not exceed 20,000. The
integers correspond to the number of couples in a population possessing each
genotype pairing for a given factor. In order, the six given integers represent
the number of couples having the following genotypes:

AA-AA
AA-Aa
AA-aa
Aa-Aa
Aa-aa
aa-aa

Return: The expected number of offspring displaying the dominant phenotype in
the next generation, under the assumption that every couple has exactly two
offspring.

Sample Dataset
--------------

1 0 0 1 0 1

Sample Output
-------------

3.5

"""
import sys

with open(sys.argv[1], 'r') as in_file:   # to automagically close the file
  p = list(map(int,in_file.readline().split()))

r = 2*(p[0]+p[1]+p[2]) + 1.5*p[3] + p[4]  # each pair has two children

print(r)
