#!/usr/bin/env python3

"""
Newick Format with Edge Weights
===============================

Given: A collection of n weighted trees (n≤40) in Newick format, with each tree
containing at most 200 nodes; each tree Tk is followed by a pair of nodes xk
and yk in Tk.

Return: A collection of n numbers, for which the kth number represents the
distance between xk and yk in Tk.

Sample Dataset
--------------

(dog:42,cat:33);
cat dog

((dog:4,cat:3):74,robot:98,elephant:58);
dog elephant

Sample Output
-------------

75 136

"""

import sys, re

def distance_in_weighted_tree(tree,a,b):
  """
  Needs a tree in the Newick format and two nodes, returns the distance
  between the nodes. It does not build the tree or do anything fancy,
  we just count parentheses representing parent nodes and find the most
  recent common ancestor, then sum both distances from a and b to the mrca.
  """
  # return zero if we got two identical nodes
  if a == b:
    return 0
  # create lists for parents
  a_parents = []
  a_distances = [0]
  b_parents = []
  b_distances = [0]
  regex = re.compile(("[0-9]+"))
  # now we search for all parents of node a, then node b
  for node, parents, distances in ([a, a_parents, a_distances], [b, b_parents, b_distances]):
    loc = tree.find(node)
    # first we check whether the node is an internal named node
    if loc > 0 and tree[loc-1] == ")":
      # if yes, we store it as the parent number zero,
      # because it might be a parent of the other node
      parents.append(loc-1)
    else:
      # if it is not an intrnal node, we put a negative number 
      parents.append(-1)
    # openings is used to store the number of "active" left parentheses
    # so we can ignore the same number of right parentheses
    # these belong to branches, not parents
    dist = regex.search(tree, loc)
    if dist:
      distances.append(dist.group())
    else:
      print("Error")
    openings = 0
    # now scan the tree on right from the node
    for i in range(loc, len(tree)):
      if tree[i] == "(":
        openings += 1
      if tree[i] == ")":
        if openings > 0:
          # ignore, it is paired with a left parenthesis
          openings -= 1
        else:
          # hey, this is a parent, store
          parents.append(i)
          dist = regex.search(tree, i)
          if dist:
            #print(dist.group())
            distances.append(dist.group())
          else:
            None
            #print(None)
  a_distances = list(map(int,a_distances))
  b_distances = list(map(int,b_distances))
  #print(a_distances,b_distances)
  # now compare both parents lists and find the most recent common ancestor
  # and return the sum of distances to it - the distance to the mrca is equal
  # to the position of mrca in a parents list, because the position 0 is
  # reserved for the starting node itself (in case it is an internal node) 
  for i in range(len(a_parents)):
    if a_parents[i] >= 0: # skip the first element if it is not internal
      if a_parents[i] in b_parents:
        #return i+b_parents.index(a_parents[i])
        #print("lol",a_distances[:i+1],b_distances[:b_parents.index(a_parents[i])+1])
        return sum(a_distances[:i+1])+sum(b_distances[:b_parents.index(a_parents[i])+1])
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

# compute the distance for each tree
distances = []
for tree, pair in zip(trees, node_pairs):
  distances.append(distance_in_weighted_tree(tree, pair[0], pair[1]))

# print the results
print(" ".join(map(str,distances)))
