#!/usr/bin/env python3

"""
Finding a Protein Motif
=======================

To allow for the presence of its varying forms, a protein motif is represented
by a shorthand as follows: [XY] means "either X or Y" and {X} means "any amino
acid except X." For example, the N-glycosylation motif is written as
N{P}[ST]{P}.

Given: At most 15 UniProt Protein Database access IDs.

Return: For each protein possessing the N-glycosylation motif, output its given
access ID followed by a list of locations in the protein string where the motif
can be found.

Sample Dataset
--------------

A2Z669
B5ZC00
P07204_TRBM_HUMAN
P20840_SAG1_YEAST

Sample Output
-------------

B5ZC00
85 118 142 306 395
P07204_TRBM_HUMAN
47 115 116 382 409
P20840_SAG1_YEAST
79 109 135 248 306 348 364 402 485 501 614

"""
import sys
from urllib.request import urlopen

with open(sys.argv[1], 'r') as in_file:
  pids = in_file.read().split()

# initialize the protein "database"
protein_db = []
# protein_db is a list of lists containing: pid, purl, pseq,
# and a list of known motif locations within the protein

for pid in pids:
  purl = "http://www.uniprot.org/uniprot/" + pid + ".fasta"
  with urlopen(purl) as fasta:
    pseq = "".join(fasta.read().decode().split("\n")[1:])
  protein_db.append([pid, purl, pseq, []])

print(protein_db)
