#!/usr/bin/env python3

"""
Catalan Numbers and RNA Secondary Structures
============================================

Given: An RNA string s having the same number of occurrences of 'A' as 'U' and
the same number of occurrences of 'C' as 'G'. The length of the string is at
most 300 bp.

Return: The total number of noncrossing perfect matchings of basepair edges in
the bonding graph of s, modulo 1,000,000.

Sample Dataset
--------------

>Rosalind_57
AUAU

Sample Output
-------------

2

"""

import sys

def count_noncrossing_perfect_matchings(rna):
  """
  Computes recursively the number of possible noncrossing perfect matchings
  in the RNA molecule. Assumes a valid RNA fulfilling the perfect matching
  criteria (should be tested before calling this function). Uses a global
  dictionary called 'rnadict' for memoization. This dictionary must exist
  and be empty before the primitive call to this function.
  """
  # if only two bases remain, and because we know the perfect matching
  # criteria are still honored, there is exactly one way to match
  # if the string is empty, we return 1 as a neutral element for
  # multiplication with the other subinterval
  l = len(rna)
  if l <= 2:
    return 1

  # do not compute the same result more than once
  if rna in rnadict:
    return rnadict[rna]

  # create a list of all positions that, when bonded with the base at the
  # first position (rna[0]), fulfill the perfect matching criteria (i.e.
  # new intervals contain the same number of 'A' as 'U' and 'C' as 'G') 
  #
  # because we always start at the first position:
  # * we do not need to care about the circular shape of the graph
  # * the same bondings are not counted repeatedly deeper in the recursion
  #   (because the first position is not included in any new subintervals)
  #
  # because we always start with a valid interval, only one of the new
  # subintervals needs to be tested
  valid_bondings = [p for p in range(1,l,2) if rna[p]==pairing[rna[0]] and
                   test_interval(rna[1:p])]

  # now count the perfect matchings - we multiply numbers for two matching
  # subintervals (as they can be freely combined) and sum numbers for
  # all subinterval pairs (i.e. valid_bondings) 
  noncrossing_perfect_matchings = 0
  for p in valid_bondings:
    matches1 = count_noncrossing_perfect_matchings(rna[1:p])
    matches2 = count_noncrossing_perfect_matchings(rna[p+1:])
    noncrossing_perfect_matchings += matches1 * matches2

  # memoize and return the number of noncrossing perfect matchings
  rnadict[rna] = noncrossing_perfect_matchings
  return noncrossing_perfect_matchings
 
def test_interval(rna):
  """
  Returns True, if the interval contains the same number of occurrences of
  'A' as 'U' and the same number of occurrences of 'C' as 'G' and the perfect
  matching is therefore possible. Returns False otherwise.
  """
  if (rna.count("A")==rna.count("U") and rna.count("C")==rna.count("G")):
    return True
  else:
    return False
  
# open a fasta file and get the sequence (uppercased, just in case)
with open(sys.argv[1], 'r') as in_file:
  lines = in_file.read().upper().splitlines()
rna = "".join(lines[1:])

# check, whether the rna fulfills the perfect matching criteria
if not test_interval(rna):
  print(0) # no perfect matchings are possible
  sys.exit() # quit correctly, with exit code 0 (success)

# define base pairing rules
pairing = {"A":"U", "U":"A", "C":"G", "G":"C"}

# initialize a dictionary for memoization
rnadict = {}

# compute and print the result
print(count_noncrossing_perfect_matchings(rna))
