#!/usr/bin/env python3
"""
GenBank Introduction
====================

Given: A genus name, followed by two dates in YYYY/M/D format.

Return: The number of Nucleotide GenBank entries for the given
genus that were published between the dates specified.

Sample Dataset
--------------

Anthoxanthum
2003/7/25
2005/12/27

Sample Output
-------------

6

"""

import sys
from Bio import Entrez

# extract data from the input file
with open(sys.argv[1], 'r') as in_file:
  genus    = in_file.readline().strip()
  date_beg = in_file.readline().strip()
  date_end = in_file.readline().strip()

# specify your email address to stop warnings
Entrez.email = "kub4@mailinator.com"

# build the search query, for more info see:
# http://www.ncbi.nlm.nih.gov/books/NBK3837/
# http://www.ncbi.nlm.nih.gov/books/NBK49540/
query = "{}[ORGN] AND {}:{}[PDAT]".format(genus, date_beg, date_end)

# create the database handle
handle = Entrez.esearch(db="nucleotide", term=query)

print(handle)
