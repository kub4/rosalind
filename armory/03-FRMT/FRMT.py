#!/usr/bin/env python3
"""
Data Formats
============

Given: A collection of n (n≤10) GenBank entry IDs.

Return: The shortest of the strings associated with the IDs in FASTA format.

Sample Dataset
--------------

FJ817486 JX069768 JX469983

Sample Output
-------------

>gi|408690371|gb|JX469983.1| Zea mays subsp. mays clone UT3343 G2-like transcription factor mRNA, partial cds
ATGATGTATCATGCGAAGAATTTTTCTGTGCCCTTTGCTCCGCAGAGGGCACAGGATAATGAGCATGCAA
GTAATATTGGAGGTATTGGTGGACCCAACATAAGCAACCCTGCTAATCCTGTAGGAAGTGGGAAACAACG
GCTACGGTGGACATCGGATCTTCATAATCGCTTTGTGGATGCCATCGCCCAGCTTGGTGGACCAGACAGA
GCTACACCTAAAGGGGTTCTCACTGTGATGGGTGTACCAGGGATCACAATTTATCATGTGAAGAGCCATC
TGCAGAAGTATCGCCTTGCAAAGTATATACCCGACTCTCCTGCTGAAGGTTCCAAGGACGAAAAGAAAGA
TTCGAGTGATTCCCTCTCGAACACGGATTCGGCACCAGGATTGCAAATCAATGAGGCACTAAAGATGCAA
ATGGAGGTTCAGAAGCGACTACATGAGCAACTCGAGGTTCAAAGACAACTGCAACTAAGAATTGAAGCAC
AAGGAAGATACTTGCAGATGATCATTGAGGAGCAACAAAAGCTTGGTGGATCAATTAAGGCTTCTGAGGA
TCAGAAGCTTTCTGATTCACCTCCAAGCTTAGATGACTACCCAGAGAGCATGCAACCTTCTCCCAAGAAA
CCAAGGATAGACGCATTATCACCAGATTCAGAGCGCGATACAACACAACCTGAATTCGAATCCCATTTGA
TCGGTCCGTGGGATCACGGCATTGCATTCCCAGTGGAGGAGTTCAAAGCAGGCCCTGCTATGAGCAAGTC
A

"""

import sys
from Bio import Entrez

# extract data from the input file
with open(sys.argv[1], 'r') as in_file:
  genbank_ids = in_file.read().split()

# specify your email address to stop warnings
Entrez.email = "kub4@mailinator.com"

# create the database handle, for parameter explanation, see:
# http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
# db = "nucleotide" or "nuccore"
# id = genbank_ids (accepts a list)
# rettype = "fasta" (to get fasta)
handle = Entrez.efetch(db="nuccore", id=genbank_ids, rettype="fasta")

# read the handle, split the lines
dbdata = handle.read().splitlines()

# create a list of fastas (each fasta as a list of lines)
fastas = []
for line in dbdata:
  if line: # ignore empty strings
    if line[0]==">":
      fastas.append([line])
    else:
      fastas[-1].append(line)

# compute the length for each fasta (the dna sequences only)
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# NOTE: from the problem description, it is not very clear
# whether we should compare complete fastas including the
# headers, or dna sequences only, this code compares dna
# sequences only because it makes more biological sense
# (the change to compare complete fastas is trivial :)
fasta_lengths = [len("".join(f[1:])) for f in fastas]

# get the index of the shortest fasta
shortest = fasta_lengths.index(min(fasta_lengths))

# print out the shortest fasta
for line in fastas[shortest]:
  print(line)
