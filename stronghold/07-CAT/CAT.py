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

"""
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
               
def make_stat(possible):
  length = len(rna) #valid for the first round only!
  lenarray = [0]*((length+2)//2)
  for p in possible:
    a, b = p
    x=min(a,b)
    y=max(a,b)
    l = y-x-1
    if l>length//2:
      l = length-l
    lenarray[l] +=1
  print(lenarray)
  
def is_in(a, interval):
  x = min(interval)
  y = max(interval)
  if a>x and a<y:
    return True
  else:
    return False

def remove_crosslinks(possible,known):
  crosslinks = []
  for p in possible:
    a, b = p
    for k in known:
      c,d=k
      if is_in(c,p) != is_in(d,p):
        crosslinks.append(p)
        break
      if a==c or a==d or b==c or b==d:
        crosslinks.append(p)
        break        
  print("crosslinks found ", len(crosslinks), crosslinks)
  return set(possible)-set(crosslinks)
"""      
     
      




def test_interval(i):
  """
  Returns True, if the interval contains the same number of occurrences of
  'A' as 'U' and the same number of occurrences of 'C' as 'G' and the perfect
  matching is therefore possible. Returns False otherwise.
  """
  if (i.count("A")==i.count("U") and i.count("C")==i.count("G")):
    return True
  else:
    return False
  
# open a file and get the sequence (uppercased, just in case)
with open(sys.argv[1], 'r') as in_file:
  lines = in_file.read().upper().splitlines()
rna = "".join(lines[1:])

# check, whether the rna fulfills the perfect matching critera
if test_interval(rna) == False:
  print("This RNA is not suitable for perfect matching!")
  print(0) # no perfect matchings are possible
  quit()

possible = find_possible_pairs(rna)
for i in range(3):
  known = known | select_known_pairs(possible)
  possible = remove_pairs(possible, known)
  possible = remove_crosslinks(possible, known)
  #print(possible)
  #print("possible", len(possible))
  print(sorted(known))
  print("known", len(known))
make_stat(possible)

