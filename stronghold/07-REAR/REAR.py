#!/usr/bin/env python3

"""
Reversal Distance
=================

A reversal of a permutation creates a new permutation by inverting some
interval of the permutation; (5,2,3,1,4), (5,3,4,1,2), and (4,1,2,3,5) are
all reversals of (5,3,2,1,4). The reversal distance between two permutations
π and σ, written drev(π,σ), is the minimum number of reversals required to
transform π into σ (this assumes that π and σ have the same length).

Given: A collection of at most 5 pairs of permutations, all of which have
length 10.

Return: The reversal distance between each permutation pair.

Sample Dataset
--------------

1 2 3 4 5 6 7 8 9 10
3 1 5 2 7 4 9 6 10 8

3 10 8 2 5 4 7 1 6 9
5 2 3 1 7 4 10 8 6 9

8 6 7 9 4 1 3 10 2 5
8 2 7 6 9 1 5 3 10 4

3 9 10 4 1 8 6 7 5 2
2 9 8 5 1 7 3 4 6 10

1 2 3 4 5 6 7 8 9 10
1 2 3 4 5 6 7 8 9 10

Sample Output
-------------

9 4 5 7 0

"""

import sys

def get_breaks(p):
  """
  Returns a tuple of "breaks" in the permutation, or in other words,
  the places where the sequence is interrupted (the neigbours do not
  differ exactly by one).
  """
  n = len(p)
  breaks = []
  for i in range(1,n):
    if abs(p[i]-p[i-1]) > 1:
      breaks.append(i)
  return tuple(breaks)

def reversal_distance(reduced):
  """
  Returns the reversal distance between 'reduced' and ordered list [1 to n],
  where 'n' is the length of 'reduced' and 'reduced' is a permutation of
  positive integers from 1 to n (inclusive), in arbitrary order.
  """
  n = len(reduced)
  assert n > 0 and set(reduced) == set(range(1,n+1)), "invalid reduced pair!"
  
  # extend 'reduced' by zero on the left and n+1 on the right,
  # so we can detect problems even on its ends
  extended = tuple([0] + reduced + [n+1])

  # now find all breaks, in other words, all the places where
  # the permutation sequence jumps by more than 1
  ext_breaks = get_breaks(extended)

  # if we have no breaks, the reversal distance is zero
  if len(ext_breaks) == 0:
    return 0

  # create a list of generations, where 'r'-th generation is a set of hopeful
  # permutations created from 'extended' by 'r' reversals, while the 0-th
  # generation contains the 'extended' permutation itself), permutations are
  # stored as tuples in format (permutation, breaks, number of breaks)
  generations = []
  first_tuple = extended,ext_breaks,len(ext_breaks)
  first_gen = set([first_tuple])
  generations.append(first_gen)
  # print(generations)
  # create the generation counter
  g = 0

  while True:
    # create statistics for the last generation
    gen_size = len(generations[g])
    gen_breaknums = [0]*(len(ext_breaks)+1)
    print("\n GENERATION", g)
    #print(gen_breaknums)
    #print(generations[g])
    for sequence, breaks, breaknum in generations[g]:
      gen_breaknums[breaknum] += 1
    print(gen_size)
    print(gen_breaknums)
    for i,x in enumerate(gen_breaknums):
      if x > 0:
        minbreak = i
        break
    # if we got the master sequence, return
    if gen_breaknums[0] > 0:
      return g
    # if we need to continue, create next generation
    next_gen = []
    for sequence, breaks, breaknum in generations[g]:
      if breaknum == minbreak:
        for i in range(0,breaknum):
          for j in range(i+1, breaknum):
            rev = sequence[:breaks[i]]+sequence[breaks[j]-1:breaks[i]-1:-1]+sequence[breaks[j]:]
            assert len(rev)==len(sequence)
            assert set(rev)==set(sequence)
            next_gen_tuple = rev, get_breaks(rev), len(get_breaks(rev))
            next_gen.append(next_gen_tuple)
            #print("GGGGGG", next_gen)
    generations.append(set(next_gen))
    # increase the generation counter
    g += 1




# read permutations from the input file
permutations = []
with open(sys.argv[1], 'r') as in_file:
  for line in in_file:
    if line.strip(): # skip empty lines
      permutations.append(line.split())

# create permutation pairs
pairs = []
assert len(permutations) % 2 == 0, "input error: uneven number of permutations"
for i,perm in enumerate(permutations):
  if i%2 == 0:
    pairs.append([perm])
  else:
    pairs[-1].append(perm)

# reduce each pair to one permutation containing integers 1...n
# (where 'n' is the length of each permutation in the pair)
# which will be compared to the ordered list [1,...,n] 
reduced_pairs = []
for p,q in pairs:
  assert len(p) == len(q),         "input error: pair length difference"
  assert len(p) == len(set(p)),    "input error: duplicate elements"
  assert set(p) == set(q),         "input error: pair element difference"
  reduced = [0]*len(p)
  for i,v in enumerate(p):
    reduced[q.index(v)] = i+1
  reduced_pairs.append(reduced)

# compute the reversal distance for each reduced pair
distances = []
for reduced in reduced_pairs:
  distances.append(reversal_distance(reduced))

# print the results
print(" ".join([str(d) for d in distances]))
