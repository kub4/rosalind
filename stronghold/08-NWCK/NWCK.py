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

import sys

def distance_in_tree(tree,a,b):
  a_parents = []
  b_parents = []
  for node, parents in ([a, a_parents], [b, b_parents]):
    loc = tree.find(node)
    if loc > 0 and tree[loc-1] == ")":
      parents.append(loc-1)
    else:
      parents.append(-1)
    openings = 0
    for i in range(loc, len(tree)):
      if tree[i] == "(":
        openings += 1
      if tree[i] == ")":
        if openings > 0:
          openings -= 1
        else:
          parents.append(i)
  # now compare both parents lists and find the most recent common ancestor
  # return the sum of distances from nodes to the most recent common ancestor
  for i in range(len(a_parents)):
    if a_parents[i] >= 0:
      if a_parents[i] in b_parents:
        return i+b_parents.index(a_parents[i])
  return None # should not happen, probably some node missing in the tree


# open the file
with open(sys.argv[1], 'r') as in_file:
  lines = in_file.readlines()

# parse the data
trees = []
node_pairs = []
for line in lines:
  line = line.strip()
  if line:
    if "(" in line: #we assume each tree contains parentheses
      trees.append(line)
    else: # but no node name contains parentheses
      node_pairs.append(line.split())

# sanity check
assert len(trees)==len(node_pairs)

# compute the distances
distances = []
for tree, pair in zip(trees, node_pairs):
  distances.append(distance_in_tree(tree, pair[0], pair[1]))

# print the results
print(" ".join(map(str,distances)))
