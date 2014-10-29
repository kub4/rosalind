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

dnacode = {
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

"""
function definitions:
"""

# from a _DNA_ string, return all ORFs as a list of protein strings
#  * does not analyze the reverse complement
#  * does not remove duplicate proteins
#  * for each found START, translation is attempted by a helper function
def find_orfs(dna):
  i = 0
  proteins = []
  while i < (len(dna)-2):
    start = dna.find("ATG", i)
    if start < 0:
      break
    else:
      prot = dna_orf_translate(dna[start:])
      if prot:
        proteins.append(prot)
    i = start + 1
  return proteins

# a _DNA_ translation function with special properties:
#  * starts from the beginning of a DNA string
#  * translates until STOP is encountered
#  * if STOP is encountered, returns a protein string
#  * of no STOP is encountered, no protein is returned!
def dna_orf_translate(dna):
  prot = []
  for i in range(0, len(dna)-2, 3):
    aa = dnacode.get(dna[i:i+3])
    if aa:
      prot.append(aa)
    else:
      return "".join(prot)
  return None


"""
end of function definitions
"""

# open the file, extract the dna string
with open(sys.argv[1], 'r') as in_file:
  dna = "".join(in_file.read().upper().split()[1:])

# create the reverse complement
basepairs = str.maketrans("ATCG","TAGC")
rev_dna   = dna[::-1].translate(basepairs)

# create a list of ORF proteins from both strands
proteins = find_orfs(dna) + find_orfs(rev_dna)

# print the results, distinct candidates only
for p in set(proteins):
  print(p)
