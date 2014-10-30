#!/usr/bin/env python3

"""
Completing a Tree
=================

An undirected graph is connected if there is a path connecting any two nodes.
A tree is a connected (undirected) graph containing no cycles; this definition
forces the tree to have a branching structure organized around a central core
of nodes, just like its living counterpart.

In the creation of a phylogeny, taxa are encoded by the tree's leaves, or nodes
having degree 1. A node of a tree having degree larger than 1 is called an
internal node.

Given: A positive integer n (n≤1000) and an adjacency list corresponding to
a graph on n nodes that contains no cycles.

Return: The minimum number of edges that can be added to the graph to produce
a tree.

Sample Dataset
--------------

10
1 2
2 8
4 10
5 9
6 10
7 9

Sample Output
-------------

3

"""
import sys

# opening a file and extracting useful data
with open(sys.argv[1], 'r') as in_file:
  data = in_file.read().strip().split("\n")

n   = int(data.pop(0))
adj = [set(map(int,i.split())) for i in data]

# now we need to group the connected nodes
groups = [set(adj[0])] #initialize

for i, a in enumerate(adj):
  for j, g in enumerate(groups):
    if a & g:
      groups[j] = a|g
      break
  else:
    groups.append(a)

# count the created groups
g_count = len(groups)

# count all grouped nodes
g_nodes = sum(len(g) for g in groups)

# now get the number of ungrouped nodes
u_count = n - g_nodes

# now all the groups, including the ungrouped nodes
tot_grp = g_count + u_count

# we need one less edge to connect all groups together
print(tot_grp-1)
