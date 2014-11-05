#!/usr/bin/env python3

"""
Longest Increasing Subsequence
==============================

A subsequence of a permutation is a collection of elements of the permutation
in the order that they appear. For example, (5, 3, 4) is a subsequence of (5,
1, 3, 4, 2). A subsequence is increasing if the elements of the subsequence
increase, and decreasing if the elements decrease.

Given: A positive integer n≤10000 followed by a permutation π of length n.

Return: A longest increasing subsequence of π, followed by
a longest decreasing subsequence of π.

Sample Dataset
--------------

5
5 1 4 2 3

Sample Output
-------------

1 2 3
5 4 2

"""

import sys
from operator import neg

# kill all the hopeless candidates
# return the list containing the longest one with the smallest last item
# (that means, with the highest hopes of further extension) + those shorter
# candidates, which have even smaller last items (progressively)
def candidate_killer(candidates):
  scoreboard = [] # here we will sort the candidates
  for i, c in enumerate(candidates): # for each candidate in the list of candidates, fill its
    scoreboard.append((len(c),neg(c[-1]),i)) # length, negative of the last item, index in the list
  scoreboard.sort(reverse=True) # reverse sort, longest with smallest last item first (because we use the negative)
  good_candidates = [candidates[scoreboard[0][2]]] # put the first one into the list of good candidates
  best_last = scoreboard[0][1] # keep the information of the smallest last item so far (we use the negative)
  for s in scoreboard: # from the scoreboard
    if s[1] > best_last: # if the last item is smaller than the best known (we compare its negative!)
      good_candidates.append(candidates[s[2]]) # append the candidate into the list of good ones
      best_last = s[1] # and update the smallest last item (negative) 
  return good_candidates # return the list of good candidates
  
# search for the longest increasing subsequence
def lgis(sequence):
  candidates = [] # here we keep candidates
  for s in sequence: # for each number in the sequence
    new_candidates = [] # create an empty list for new candidates
    for c in candidates: # for each existing candidate
      if s > c[-1]: # try to create a new candidate by extension
        new_candidates.append(c+[s]) # (keeping the old copy untouched)
    new_candidates.append([s]) # use the actual number as a new, one-item candidate
    candidates += new_candidates # combine old and new candidates into one list
    candidates = candidate_killer(candidates) # keep only good candidates, kill the rest
  return max(candidates,key=len) # return the longest surviving candidate

# open a file and parse the data
with open(sys.argv[1], 'r') as in_file:
  lines = in_file.read().splitlines()

n = int(lines[0].strip())
sequence = list(map(int,lines[1].split()))
assert n == len(sequence), "Invalid input!"

# start search; for decreasing sequence, we use the same code,
# but negate everything
longest_increasing = lgis(sequence)
longest_decreasing = map(neg,lgis(map(neg,sequence)))

print(" ".join(map(str,longest_increasing)))
print(" ".join(map(str,longest_decreasing)))
