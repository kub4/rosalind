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

def find_possible_pairs(rna):
  possible_pairs = []
  l = len(rna)
  for base in [("A","U"),("C","G")]:
    for i, nt in enumerate(rna):
      if nt == base[0]:
        for j, nu in enumerate(rna):
          if nu == base[1]:
            x = min(i,j)
            y = max(i,j)
            frg1 = rna[x+1:y]
            frg2 = rna[y+1:]+rna[:x]
            if test(frg1) and test(frg2):
              possible_pairs.append((i,j))
  return possible_pairs

def select_known_pairs(possible):
  # count how many times occur each rna position in pairs
  counts = [0]*len(rna)
  for p in possible:
    counts[p[0]] += 1
    counts[p[1]] += 1
  print(counts)
  # select those occuring only once (we hope for them!)
  magic_pos = []
  for i, val in enumerate(counts):
    if val == 1:
      magic_pos.append(i)
  #print(magic_pos)
  # find the known pairs
  known_pairs = []
  for m in magic_pos:
    for p in possible:
      if m == p[0] or m == p[1]:
        known_pairs.append(p)
  return set(known_pairs)
               

def remove_pairs(masterlist,remove):
  unique = set(masterlist)-set(remove)
  return unique  

def test(interval):
  a=len([a for a in interval if a == "A"])
  c=len([c for c in interval if c == "C"])
  g=len([g for g in interval if g == "G"])
  u=len([u for u in interval if u == "U"])
  if a==u and c==g:
    return True
  else:
    return False
  

# open a file and get the sequence
with open(sys.argv[1], 'r') as in_file:
  lines = in_file.read().upper().splitlines()
rna = "".join(lines[1:])

possible = set()
known = set()

possible = find_possible_pairs(rna)
for i in range(10):
  known = known | select_known_pairs(possible)
  possible = remove_pairs(possible, known)
  #print(possible)
  #print("possible", len(possible))
  #print(known)
  print("known", len(known))

