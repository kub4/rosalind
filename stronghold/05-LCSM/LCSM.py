#!/usr/bin/env python3

"""
Finding a Shared Motif
======================

A common substring of a collection of strings is a substring of every member of
the collection. We say that a common substring is a longest common substring if
there does not exist a longer common substring. For example, "CG" is a common
substring of "ACGTACGT" and "AACCGGTATA", but it is not as long as possible; in
this case, "GTA" is a longest common substring of "ACGTACGT" and "AACCGTATA".

The longest common substring is not necessarily unique; for a simple example,
"AA" and "CC" are both longest common substrings of "AACC" and "CCAA".

Given: A collection of k (k≤100) DNA strings of length at most 1 kbp each in
FASTA format.

Return: A longest common substring of the collection. (If multiple solutions
exist, you may return any single solution.)

Sample Dataset
--------------

>Rosalind_1
GATTACA
>Rosalind_2
TAGACCA
>Rosalind_3
ATACA

Sample Output
-------------

AC

"""
import sys

# define a function generating all substrings
# starting with the longest (the string itself)
def subs(string):
  l = len(string)
  for i in range(l):
    for j in range(0,i+1):
      yield string[j:l-i+j]
"""
end of function definitions
"""

with open(sys.argv[1], 'r') as in_file:   # to automagically close the file
  data = in_file.read().strip()

raw_seq = filter(None,data.split(">"))    # filter removes empty strings
sequences = []                            # a list of lists [[sequence, length]]

for s in raw_seq:                         # parsing the sequences
  s = s.split()
  sdata = ''.join(s[1:])
  sequences.append([sdata, len(sdata)])

"""
ok, now we have a nice list of fasta sequences,
including lengths, now to search through the haystack...
"""

# choose one shortest sequence and create a set of other sequences
shortseq = min(sequences, key = lambda x: x[1])[0]
otherseq = set([x[0] for x in sequences if x[0] != shortseq])
"""
now we will create all subsequences of the shortseq (starting from
the longest, the shortseq itself) and test their presence in other
sequences...
"""

for sub in subs(shortseq):
  lcs = True # is it the longest common sequence? initialized as True
  for other in otherseq:
    if other.find(sub)==-1:
      lcs = False # if not found, change lcs to False
      break
  if lcs==1:
    print(sub) # if lcs still holds True, print the result
    break      # and end the search (we need only one)
