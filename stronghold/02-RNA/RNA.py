#!/usr/bin/env python3

"""
Transcribing DNA into RNA
=========================

Given: A DNA string t corresponding to a coding strand (at most 1000 nt).

Return: The transcribed RNA string of t.

Sample Dataset
--------------

GATGGAACTTGACTACGTAAATT

Sample Output
-------------

GAUGGAACUUGACUACGUAAAUU

"""

import sys

with open(sys.argv[1], 'r') as input_file:   # to automagically close the file
  data = input_file.read()                   # when leaving the nested block

dna = ''.join(data.split())  # the split/join trick to remove all whitespace

transcription_table = str.maketrans("Tt","Uu")  # preserving case
rna = dna.translate(transcription_table)        # much confusement, very chaos!

print(rna)
