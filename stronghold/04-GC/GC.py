#!/usr/bin/env python3

"""
Computing GC Content
====================
Given: At most 10 DNA strings in FASTA format (of length at most 1 kbp each).

Return: The ID of the string having the highest GC-content, followed by the
GC-content of that string. Rosalind allows for a default error of 0.001 in all
decimal answers unless otherwise stated.

Sample Dataset
--------------

>Rosalind_6404
CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC
TCCCACTAATAATTCTGAGG
>Rosalind_5959
CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCT
ATATCCATTTGTCAGCAGACACGC
>Rosalind_0808
CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGAC
TGGGAACCTGCGGGCAGTAGGTGGAAT

Sample Output
-------------

Rosalind_0808
60.919540

"""
import sys

def gccontent(dna):
  dna = dna.upper()
  return (dna.count("G") + dna.count("C"))/len(dna)

with open(sys.argv[1], 'r') as in_file:   # to automagically close the file
  data = in_file.read().strip()

raw_seq = filter(None,data.split(">"))    # filter removes empty strings
sequences = []        # a list of lists [[name, gccontent, sequence]]

for s in raw_seq:
  s = s.split()
  sequences.append([s[0],0,''.join(s[1:])])

maxgc = 0
for s in sequences:
  s[1] = gccontent(s[2])
  if s[1]>maxgc: maxgc=s[1]

for s in sequences: # in case of a tie, print all
  if s[1]==maxgc:
    print(s[0])
    print(s[1]*100)
