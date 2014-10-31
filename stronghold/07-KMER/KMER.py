#!/usr/bin/env python3

"""
k-Mer Composition
=================

For a fixed positive integer k, order all possible k-mers taken from
an underlying alphabet lexicographically. Then the k-mer composition
of a string s can be represented by an array A for which A[m] denotes
the number of times that the mth k-mer (with respect to the lexicographic
order) appears in s.

Given: A DNA string s in FASTA format (having length at most 100 kbp).

Return: The 4-mer composition of s.


Sample Dataset
--------------

>Rosalind_6431
CTTCGAAAGTTTGGGCCGAGTCTTACAGTCGGTCTTGAAGCAAAGTAACGAACTCCACGG
CCCTGACTACCGAACCAGTTGTGAGTACTCAACTGGGTGAGAGTGCAGTCCCTATTGAGT
TTCCGAGACTCACCGGGATTTTCGATCCAGCCTCAGTCCAGTCTTGTGGCCAACTCACCA
AATGACGTTGGAATATCCCTGTCTAGCTCACGCAGTACTTAGTAAGAGGTCGCTGCAGCG
GGGCAAGGAGATCGGAAAATGTGCTCTATATGCGACTAAAGCTCCTAACTTACACGTAGA
CTTGCCCGTGTTAAAAACTCGGCTCACATGCTGTCTGCGGCTGGCTGTATACAGTATCTA
CCTAATACCCTTCAGTTCGCCGCACAAAAGCTGGGAGTTACCGCGGAAATCACAG

Sample Output
-------------

4 1 4 3 0 1 1 5 1 3 1 2 2 1 2 0 1 1 3 1 2 1 3 1 1 1 1 2 2 5 1 3 0 2 2 1 1 1 1 3 1 0 0 1 5 5 1 5 0 2 0 2 1 2 1 1 1 2 0 1 0 0 1 1 3 2 1 0 3 2 3 0 0 2 0 8 0 0 1 0 2 1 3 0 0 0 1 4 3 2 1 1 3 1 2 1 3 1 2 1 2 1 1 1 2 3 2 1 1 0 1 1 3 2 1 2 6 2 1 1 1 2 3 3 3 2 3 0 3 2 1 1 0 0 1 4 3 0 1 5 0 2 0 1 2 1 3 0 1 2 2 1 1 0 3 0 0 4 5 0 3 0 2 1 1 3 0 3 2 2 1 1 0 2 1 0 2 2 1 2 0 2 2 5 2 2 1 1 2 1 2 2 2 2 1 1 3 4 0 2 1 1 0 1 2 2 1 1 1 5 2 0 3 2 1 1 2 2 3 0 3 0 1 3 1 2 3 0 2 1 2 2 1 2 3 0 1 2 3 1 1 3 1 0 1 1 3 0 2 1 2 2 0 2 1 1

"""

import sys
from itertools import product
from collections import Counter

# open a file and get a dna sequence
with open(sys.argv[1], 'r') as in_file:
  data = in_file.read().split()

dna = ''.join(data[1:]).upper()

# we want 4-mers
k = 4

# count the k-mers in the dna string
dna_kmers = Counter([dna[i:i+k] for i in range(len(dna)-k+1)])

# generate all possible k-mers in lexicographic order
possible_kmers = ["".join(x) for x in product("ACGT", repeat=k)]

# build the lexicographical count list for all possible k-mers
count_list = [dna_kmers[kmer] for kmer in possible_kmers]

# format and print
print(" ".join(map(str, count_list)))
