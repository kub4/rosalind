#!/usr/bin/env python3

"""
RNA Splicing
============

Given: A DNA string s (of length at most 1 kbp) and a collection of substrings
of s acting as introns. All strings are given in FASTA format.

Return: A protein string resulting from transcribing and translating the exons
of s. (Note: Only one solution will exist for the dataset provided.)

Sample Dataset
--------------

>Rosalind_10
ATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTCGAGAGGCATATGGTCACATGATCGGTCGAGCGTGTTTCAAAGTTTGCGCCTAG
>Rosalind_12
ATCGGTCGAA
>Rosalind_15
ATCGGTCGAGCGTGT

Sample Output
-------------

MVYIADKQHVASREAYGHMFKVCA

"""
import sys

code = {
"UUU": "F",    "CUU": "L",    "AUU": "I",    "GUU": "V",
"UUC": "F",    "CUC": "L",    "AUC": "I",    "GUC": "V",
"UUA": "L",    "CUA": "L",    "AUA": "I",    "GUA": "V",
"UUG": "L",    "CUG": "L",    "AUG": "M",    "GUG": "V",
"UCU": "S",    "CCU": "P",    "ACU": "T",    "GCU": "A",
"UCC": "S",    "CCC": "P",    "ACC": "T",    "GCC": "A",
"UCA": "S",    "CCA": "P",    "ACA": "T",    "GCA": "A",
"UCG": "S",    "CCG": "P",    "ACG": "T",    "GCG": "A",
"UAU": "Y",    "CAU": "H",    "AAU": "N",    "GAU": "D",
"UAC": "Y",    "CAC": "H",    "AAC": "N",    "GAC": "D",
"UAA": "STOP", "CAA": "Q",    "AAA": "K",    "GAA": "E",
"UAG": "STOP", "CAG": "Q",    "AAG": "K",    "GAG": "E",
"UGU": "C",    "CGU": "R",    "AGU": "S",    "GGU": "G",
"UGC": "C",    "CGC": "R",    "AGC": "S",    "GGC": "G",
"UGA": "STOP", "CGA": "R",    "AGA": "R",    "GGA": "G",
"UGG": "W",    "CGG": "R",    "AGG": "R",    "GGG": "G"
}

#translate RNA from AUF to STOP or end of RNA
def translate(rna):
  length = len(rna)
  start  = rna.find("AUG")
  protein = []
  if start >= 0:
    for i in range(start, length-2, 3):
      aa = code.get(rna[i:i+3])
      if aa == "STOP":
        break
      protein.append(aa)
  return "".join(protein)
      
"""
end of functions definitions
"""
# read the file, extract the RNA and a list of introns
with open(sys.argv[1], 'r') as in_file:
  data = in_file.read().upper().replace("T","U").strip()
  # when splicing, better work with RNA :)

raw_seq = filter(None,data.split(">"))  # filter removes empty strings
sequences = []        # a list of sequences (no need to keep the name)
for s in raw_seq:
  s = s.split()
  sequences.append(''.join(s[1:]))

rna = sequences[0]
introns = sequences[1:]

# remove introns...
for intron in introns:
  rna = rna.replace(intron,"")

# translate and print...
print(translate(rna))
