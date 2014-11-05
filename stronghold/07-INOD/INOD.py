#!/usr/bin/env python3

"""
Counting Phylogenetic Ancestors
===============================

Even though a binary tree can include nodes having degree 2, an unrooted binary
tree is defined more specifically: all internal nodes have degree 3. In turn,
a rooted binary tree is such that only the root has degree 2 (all other
internal nodes have degree 3).

Given: A positive integer n (3≤n≤10000).

Return: The number of internal nodes of any unrooted binary tree
having n leaves.

Sample Dataset
--------------

4

Sample Output
-------------

2

"""

import sys
#from math import factorial as fac

# open a file and get n
with open(sys.argv[1], 'r') as in_file:
  n = int(in_file.read().strip())

# after drawing some trees, answer is very simple:
# number of internal nodes = leaves - 2
# and if you add another leaf, you create one more
# internal node... so...

print(n-2)
