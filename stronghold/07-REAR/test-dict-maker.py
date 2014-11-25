#!/usr/bin/env python3

"""
A simple script to create a dictionary containing all permutations
of range(1,11) (as keys) and their reversal distances from the
ordered permutation (1, 2, 3, 4, 5, 6, 7, 8, 9, 10). The created
dictionary is intended for testing Rosalind's REAR algorithms.
WARNING: the dictionary will contain 10! (3628800) keys/values!
"""
from math import factorial as fac
def get_revdist1(p):
  length = len(p)
  reversions = []
  for i in range(length-1):
    for j in range(i+1,length):
      rev = p[:i] + tuple(reversed(p[i:j+1])) + p[j+1:]
      assert len(rev) == len(p) and set(rev) == set(p)
      reversions.append(rev)
  return reversions

# create the permutations dictionary, add the initial ordered permutation
permutations = {tuple(range(1,11)):0}

# create the last generation dictionary
last_generation = dict(permutations)

# initialize the generation counter
g = 0

while True:
  print(g,len(permutations),len(permutations)/fac(10)*100,"%")
  g += 1
  new_generation = {}
  for p in last_generation:
    revdist1 = get_revdist1(p) # a list of all reversions of distance 1
    for r in revdist1:
      if r not in permutations:
        new_generation[r] = g
  if len(new_generation) == 0:
    break
  permutations.update(new_generation)
  last_generation = dict(new_generation)

#print(permutations)
