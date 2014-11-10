#!/usr/bin/env python3

"""
Error Correction in Reads
=========================

As is the case with point mutations, the most common type of sequencing error
occurs when a single nucleotide from a read is interpreted incorrectly.

Given: A collection of up to 1000 reads of equal length (at most 50 bp) in
FASTA format. Some of these reads were generated with a single-nucleotide
error. For each read s in the dataset, one of the following applies:

 * s was correctly sequenced and appears in the dataset at least twice
   (possibly as a reverse complement);

 * s is incorrect, it appears in the dataset exactly once, and its
   Hamming distance is 1 with respect to exactly one correct read
   in the dataset (or its reverse complement).

Return: A list of all corrections in the form "[old read]->[new read]".
(Each correction must be a single symbol substitution, and you may return
the corrections in any order.)

Sample Dataset
--------------

>Rosalind_52
TCATC
>Rosalind_44
TTCAT
>Rosalind_68
TCATC
>Rosalind_28
TGAAA
>Rosalind_95
GAGGA
>Rosalind_66
TTTCA
>Rosalind_33
ATCAA
>Rosalind_21
TTGAT
>Rosalind_18
TTTCC

Sample Output
-------------

TTCAT->TTGAT
GAGGA->GATGA
TTTCC->TTTCA

"""

import sys
from collections import defaultdict

def reverse_complement(dna):
  """
  Returns the reverse complement DNA string. Requires an uppercase DNA string.
  Requires 'revdict' dictionary for memoization.
  """
  if dna not in revdict:
    revdict[dna] = dna[::-1].translate(str.maketrans("ACGT","TGCA"))
  return revdict[dna]

def ham1_compare(dna1,dna2):
  """
  Compares two strings, returns True if the Hamming distance
  between them is 0 or 1. Returns False otherwise.
  """ 
  diff = 0
  for a,b in zip(dna1,dna2):
    if a != b:
      diff += 1
      if diff > 1:
        return False
  return True

# open the file and get the sequences
with open(sys.argv[1], 'r') as in_file:
  lines = in_file.readlines()

# open the input file and get the sequences, uppercased just in case
sequences = []
for line in lines:
  if line[0]==">":
    sequences.append([])
  elif sequences:
    sequences[-1].append(line.strip().upper())
sequences = ["".join(s) for s in sequences]

# create a default dictionary for counting unique sequences, treating
# reverse complements as synonyms, storing always the version that is
# the first in lexicographical order
uniques = defaultdict(int)

# creates a dictionary for reverse complement memoization
revdict = {}

# create lists to hold all bad reads and their corrected versions
bad_reads = []
corrected_reads = []

# count unique sequences, treating reverse complements as synonyms,
# storing the version that is the first in lexicographical order
for s in sequences:
  r = reverse_complement(s)
  if s < r:
    uniques[s] += 1
  else:
    uniques[r] += 1

# identify all bad reads and add them to their own list, if the bad read
# is not in the original seqeunces list, add its reverse complement instead
for uniq, num in uniques.items():
  if num == 1:
    if uniq in sequences:
      bad_reads.append(uniq)
    else:
      bad_reads.append(reverse_complement(uniq))

# find corrected versions for all bad reads, add them to the list
# compare also with reverse complements and use the appropriate version
for bad in bad_reads:
  for uniq, num in uniques.items():
    if num >= 2:
      if ham1_compare(bad,uniq):
        corrected_reads.append(uniq)
        break # we need only one correct version
      if ham1_compare(bad,reverse_complement(uniq)):
        corrected_reads.append(reverse_complement(uniq))
        break # we need only one correct version

# print out the results
for bad, corr in zip(bad_reads, corrected_reads):
  print(bad,"->",corr,sep="")
