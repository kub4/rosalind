#!/usr/bin/env python3

"""
Overlap Graphs
==============

For a collection of strings and a positive integer k, the overlap graph for the
strings is a directed graph Ok in which each string is represented by a node,
and string s is connected to string t with a directed edge when there is a
length k suffix of s that matches a length k prefix of t, as long as s≠t; we
demand s≠t to prevent directed loops in the overlap graph (although directed
cycles may be present).

Given: A collection of DNA strings in FASTA format having total length at most
10 kbp.

Return: The adjacency list corresponding to O3. You may return edges in any
order.

Sample Dataset
--------------

>Rosalind_0498
AAATAAA
>Rosalind_2391
AAATTTT
>Rosalind_2323
TTTTCCC
>Rosalind_0442
AAATCCC
>Rosalind_5013
GGGTGGG

Sample Output
-------------

Rosalind_0498 Rosalind_2391
Rosalind_0498 Rosalind_0442
Rosalind_2391 Rosalind_2323

"""
import sys


with open(sys.argv[1], 'r') as in_file:   # to automagically close the file
  data = in_file.read().strip()

raw_seq = filter(None,data.split(">"))    # filter removes empty strings
sequences = []                            # a list of lists [[name, sequence]]

for s in raw_seq:
  s = s.split()
  sequences.append([s[0],''.join(s[1:])])

"""
ok, now we have a nice list of fasta sequences,
now to analyze adjacencies...
"""

for s in sequences:
  for t in sequences:
    if s[1][-3:]==t[1][:3] and s[0]!=t[0]:
      print(s[0]+' '+t[0])
