#!/usr/bin/env python3

"""
Inferring mRNA from Protein
===========================

Given: A protein string of length at most 1000 aa.

Return: The total number of different RNA strings from which the protein could
have been translated, modulo 1,000,000. (Don't neglect the importance of the
stop codon in protein translation.)

Sample Dataset
--------------

MA

Sample Output
-------------

12

"""
import sys
from collections import Counter

code = {
"UUU" : "F", "CUU" : "L", "AUU" : "I", "GUU" : "V",
"UUC" : "F", "CUC" : "L", "AUC" : "I", "GUC" : "V",
"UUA" : "L", "CUA" : "L", "AUA" : "I", "GUA" : "V",
"UUG" : "L", "CUG" : "L", "AUG" : "M", "GUG" : "V",
"UCU" : "S", "CCU" : "P", "ACU" : "T", "GCU" : "A",
"UCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A",
"UCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A",
"UCG" : "S", "CCG" : "P", "ACG" : "T", "GCG" : "A",
"UAU" : "Y", "CAU" : "H", "AAU" : "N", "GAU" : "D",
"UAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
"UAA" : "STOP", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
"UAG" : "STOP", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
"UGU" : "C", "CGU" : "R", "AGU" : "S", "GGU" : "G",
"UGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
"UGA" : "STOP", "CGA" : "R", "AGA" : "R", "GGA" : "G",
"UGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G"
}

values = []

for c in code:
  values.append(code[c])
c = Counter(values)

for key, value in sorted(c.items()):
  print('"'+str(key)+'"'+' : '+str(value))
   
