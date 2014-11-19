#!/usr/bin/env python3
"""
Introduction to Protein Databases
=================================

Given: The UniProt ID of a protein.

Return: A list of biological processes in which the protein is involved
(biological processes are found in a subsection of the protein's "Gene
Ontology" (GO) section).

Sample Dataset
--------------

Q5SLP9

Sample Output
-------------

DNA recombination
DNA repair
DNA replication

"""
import sys
from Bio import ExPASy
from Bio import SwissProt

# extract the protein id from the input file
with open(sys.argv[1], 'r') as in_file:
  pid = in_file.read().strip()

# create a handle...
# note that the comment "you can give several IDs separated by commas"
# in the rosalind's programming shorcut instructions is misleading, tried
# it during SWAT (where two protein ids are given) and it does not work...
# and after a short look into the biopython's source, it cannot work...
handle = ExPASy.get_sprot_raw(pid)

# read it
record = SwissProt.read(handle)

"""
Now we need to parse... The GO section is a part of DR (Database
cross-References), so we can simply ask for 'cross_references'...
and select the GO lines. For GO, the third field is a 1-letter
abbreviation for one of the 3 ontology aspects (P for biological
Process, F for molecular Function and C for cellular Component),
separated from the GO term by a column... We need only biological
processes. See http://web.expasy.org/docs/userman.html#DR_line
"""
# get database cross-references (DR lines), keeping only GO lines
xrefs = [x for x in record.cross_references if x[0]=="GO"]

# select the biological processes, cutting the [PFC]: part
procs = [x[2][2:] for x in xrefs if x[2][0]=="P"]

# print the results
print(*procs, sep='\n')
