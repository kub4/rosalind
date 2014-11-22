#!/usr/bin/env python3
"""
Protein Translation
===================

The Translate tool from the SMS 2 package can be found in the SMS 2 package:
http://www.bioinformatics.org/sms2/translate.html

A detailed list of genetic code variants (codon tables) along with indexes
representing these codes (1 = standard genetic code, etc.) is here:
http://www.bioinformatics.org/JaMBW/2/3/TranslationTables.html

For now, when translating DNA and RNA strings, we will start with the first
letter of the string and ignore stop codons.

Given: A DNA string s of length at most 10 kbp, and a protein translated by s.

Return: The index of the genetic code variant that was used for translation.
(If multiple solutions exist, you may return any one.)

Sample Dataset
--------------

ATGGCCATGGCGCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA
MAMAPRTEINSTRING

Sample Output
-------------

1

"""

import sys
from Bio.Seq import translate
from Bio.Seq import CodonTable

# extract and purify biopolymers from the input file
with open(sys.argv[1], "r") as in_file:
  dna            = in_file.readline().strip()
  master_protein = in_file.readline().strip()

"""
# manually taken list from CodonTable.py (Biopython source),
# the list at http://www.bioinformatics.org/JaMBW/2/3/TranslationTables.html
# is incomplete (Last update of the Genetic Codes: Sep 26, 1996)
valid_tables = [1,2,3,4,5,6,9,10,11,12,13,14,15,16,21,22,23]
"""
# ok, it is better to get the valid tables list programatically
# print(CodonTable.unambiguous_dna_by_id)
valid_tables = [k for k,v in CodonTable.unambiguous_dna_by_id.items()]

# a list of the codes possibly used for translating our dna to our protein
# (yet empty)
used_codes = []

# now we translate using all valid tables and check whether the resulting
# protein is the same as our given master protein
for t in valid_tables:
  # 'stop_symbol=""' and 'to_stop=False' to IGNORE STOP CODONS
  # 'cds=False' to ignore coding sequence checking (whether the sequence
  # starts with START, whether the sequence length is a multiple of three...
  protein = translate(dna, table=t, stop_symbol="", to_stop=False, cds=False)
  if protein == master_protein:
    used_codes.append(t)

# if we had found some possible codes for our protein, print the first one
# otherwise, print None
if used_codes:
  print(used_codes[0])
else:
  print(None)
