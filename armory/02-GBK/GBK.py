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

# the following code would, from the handle, read a xml-formatted
# string, that would need to be parsed (note also that the handle
# is closed after use, so you cannot open it once and then compare
# two ways of reading the handle without reopening it - therefore
# simple uncommenting of the following line would result in
# "OSError: Can't parse a closed handle" comming from the second
# read attempt)

###record = handle.read()

# however the following code parses the xml string automagically,
# record then contains <class 'Bio.Entrez.Parser.DictionaryElement'>,
# which can be inspected using 'print(record') and contains, among
# other items "'Count': '6'"...

record = Entrez.read(handle)

# so we print out the value for 'Count'
print(record["Count"])
