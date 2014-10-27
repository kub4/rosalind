#!/usr/bin/env python3

"""
Translating RNA into Protein
============================

The RNA codon table dictates the details regarding the encoding of specific
codons into the amino acid alphabet.

UUU F      CUU L      AUU I      GUU V
UUC F      CUC L      AUC I      GUC V
UUA L      CUA L      AUA I      GUA V
UUG L      CUG L      AUG M      GUG V
UCU S      CCU P      ACU T      GCU A
UCC S      CCC P      ACC T      GCC A
UCA S      CCA P      ACA T      GCA A
UCG S      CCG P      ACG T      GCG A
UAU Y      CAU H      AAU N      GAU D
UAC Y      CAC H      AAC N      GAC D
UAA Stop   CAA Q      AAA K      GAA E
UAG Stop   CAG Q      AAG K      GAG E
UGU C      CGU R      AGU S      GGU G
UGC C      CGC R      AGC S      GGC G
UGA Stop   CGA R      AGA R      GGA G
UGG W      CGG R      AGG R      GGG G

Given: An RNA string s corresponding to a strand of mRNA
(of length at most 10 kbp).

Return: The protein string encoded by s.

Sample Dataset
--------------

AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA

Sample Output
-------------

MAMAPRTEINSTRING

"""
import sys

with open(sys.argv[1], 'r') as in_file:   # to automagically close the file
  mrna = in_file.read().strip()           # when leaving the nested block

code = {
"UUU" : "F",  "CUU" : "L",  "AUU" : "I",  "GUU" : "V",
"UUC" : "F",  "CUC" : "L",  "AUC" : "I",  "GUC" : "V",
"UUA" : "L",  "CUA" : "L",  "AUA" : "I",  "GUA" : "V",
"UUG" : "L",  "CUG" : "L",  "AUG" : "M",  "GUG" : "V",
"UCU" : "S",  "CCU" : "P",  "ACU" : "T",  "GCU" : "A",
"UCC" : "S",  "CCC" : "P",  "ACC" : "T",  "GCC" : "A",
"UCA" : "S",  "CCA" : "P",  "ACA" : "T",  "GCA" : "A",
"UCG" : "S",  "CCG" : "P",  "ACG" : "T",  "GCG" : "A",
"UAU" : "Y",  "CAU" : "H",  "AAU" : "N",  "GAU" : "D",
"UAC" : "Y",  "CAC" : "H",  "AAC" : "N",  "GAC" : "D",
"UAA" : None, "CAA" : "Q",  "AAA" : "K",  "GAA" : "E",
"UAG" : None, "CAG" : "Q",  "AAG" : "K",  "GAG" : "E",
"UGU" : "C",  "CGU" : "R",  "AGU" : "S",  "GGU" : "G",
"UGC" : "C",  "CGC" : "R",  "AGC" : "S",  "GGC" : "G",
"UGA" : None, "CGA" : "R",  "AGA" : "R",  "GGA" : "G",
"UGG" : "W",  "CGG" : "R",  "AGG" : "R",  "GGG" : "G"
}
