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

n   = data.pop(0)
adj = [i.split() for i in data]

# now we need to group the connected nodes
groups = [set(adj[0])] #initialize

for a in adj:
  for g in groups:
  

  

print(n, adj, groups)
