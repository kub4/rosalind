#!/usr/bin/env python3

"""
Open Reading Frames
===================

Either strand of a DNA double helix can serve as the coding strand for RNA
transcription. Hence, a given DNA string implies six total reading frames, or
ways in which the same region of DNA can be translated into amino acids: three
reading frames result from reading the string itself, whereas three more result
from reading its _reverse complement_.

Given: A DNA string s of length at most 1 kbp in FASTA format.

Return: Every distinct candidate protein string that can be translated from
ORFs of s. Strings can be returned in any order.

Sample Dataset
--------------

>Rosalind_99
AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG

Sample Output
-------------

MLLGSFRLIPKETLIQVAGSSPCNLS
M
MGMTPRLGLESLLE
MTPRLGLESLLE

"""
import sys

code = {
"TTT" : "F",  "CTT" : "L",  "ATT" : "I",  "GTT" : "V",
"TTC" : "F",  "CTC" : "L",  "ATC" : "I",  "GTC" : "V",
"TTA" : "L",  "CTA" : "L",  "ATA" : "I",  "GTA" : "V",
"TTG" : "L",  "CTG" : "L",  "ATG" : "M",  "GTG" : "V",
"TCT" : "S",  "CCT" : "P",  "ACT" : "T",  "GCT" : "A",
"TCC" : "S",  "CCC" : "P",  "ACC" : "T",  "GCC" : "A",
"TCA" : "S",  "CCA" : "P",  "ACA" : "T",  "GCA" : "A",
"TCG" : "S",  "CCG" : "P",  "ACG" : "T",  "GCG" : "A",
"TAT" : "Y",  "CAT" : "H",  "AAT" : "N",  "GAT" : "D",
"TAC" : "Y",  "CAC" : "H",  "AAC" : "N",  "GAC" : "D",
"TAA" : None, "CAA" : "Q",  "AAA" : "K",  "GAA" : "E",
"TAG" : None, "CAG" : "Q",  "AAG" : "K",  "GAG" : "E",
"TGT" : "C",  "CGT" : "R",  "AGT" : "S",  "GGT" : "G",
"TGC" : "C",  "CGC" : "R",  "AGC" : "S",  "GGC" : "G",
"TGA" : None, "CGA" : "R",  "AGA" : "R",  "GGA" : "G",
"TGG" : "W",  "CGG" : "R",  "AGG" : "R",  "GGG" : "G"
}

# the translation function, takes care of STOPs automatically
def dna_translate(dna):
  prot = []
  stop = 0
  for i in range(0, length-2, 3):
    aa = code.get(dna[i:i+3])
    print(aa)
    if aa == None:
      print("none aa")
      return "".join(prot)
    else:
      prot.append(aa)
  return None

with open(sys.argv[1], 'r') as in_file:
  dna = "".join(in_file.read().upper().split()[1:])

# measure the dna strand
length    = len(dna)

# create the reverse complement
basepairs = str.maketrans("ATCG","TAGC")
rev_dna   = dna[::-1].translate(basepairs)

# initialize the list of proteins
proteins  = []

# go through both strands and find start codons
# if found, send the relevant part of the strand
# to translation (no need to search for STOPs now)
for strand in [dna, rev_dna]:
  i = 0
  while i < (length-2):
    print(i)
    start = strand.find("ATG", i)
    if start < 0:
      print("BREAK")
      break
    else:
      print(strand[start:], start, dna_translate(strand[start:])) #DEBUG FIXME
      protein = dna_translate(strand[start:])
      if protein:
        proteins.append(protein)
    i = start + 1

print(proteins)
