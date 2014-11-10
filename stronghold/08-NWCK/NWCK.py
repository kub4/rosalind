#!/usr/bin/env python3

"""
Distances in Trees
==================

Given: A collection of n trees (n≤40) in Newick format, with each tree
containing at most 200 nodes; each tree Tk is followed by a pair of nodes xk
and yk in Tk.

Return: A collection of n positive integers, for which the kth integer
represents the distance between xk and yk in Tk.

Sample Dataset
--------------

(cat)dog;
dog cat

(dog,cat);
dog cat

Sample Output
-------------

1 2

"""

import sys, ast

def distance_in_tree(tree,a,b):
  a_parents = []
  b_parents = []
  a_loc = tree.find(a)
  b_loc = tree.find(b)
  if a_loc > 0 and tree[a_loc-1] == ")":
    a_parents.append(a_loc-1)
  else:
    a_parents.append(-100)
  if b_loc > 0 and tree[b_loc-1] == ")":
    b_parents.append(b_loc-1)
  else:
    b_parents.append(-100)
  openings = 0
  for i in range(a_loc, len(tree)):
    if tree[i] == "(":
      openings += 1
    if tree[i] == ")":
      if openings > 0:
        openings -= 1
        #print("xa",i, openings)
      else:
        #print(openings)
        a_parents.append(i)
  openings  = 0
  for i in range(b_loc, len(tree)):
    if tree[i] == "(":
      openings += 1
    if tree[i] == ")":
      if openings > 0:
        openings -= 1
        #print("xb",i, openings)
      else:
        #print(openings)
        b_parents.append(i)
  #print(a_parents)
  #print(b_parents)
  for i in range(len(a_parents)):
    if a_parents[i] >= 0:
       if a_parents[i] in b_parents:
         #print(a_parents, b_parents)
         #print(i,b_parents.index(a_parents[i]))
         return i+b_parents.index(a_parents[i])
  return 10000000


# open a file and parse the data
with open(sys.argv[1], 'r') as in_file:
  lines = in_file.readlines()

trees = []
node_pairs = []
for line in lines:
  line = line.strip()
  if line:
    if "(" in line:
      trees.append(line)
    else:
      node_pairs.append(line.split())

# sanity check
assert len(trees)==len(node_pairs)

# compute and print the distances
distances = []
for tree, pair in zip(trees, node_pairs):
  #print(tree, pair)
  distances.append(distance_in_tree(tree, pair[0], pair[1]))

print(" ".join(map(str,distances)))
