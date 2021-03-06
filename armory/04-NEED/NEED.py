#!/usr/bin/env python3
"""
Pairwise Global Alignment
=========================

Use:
 * http://www.ebi.ac.uk/Tools/psa/emboss_needle/nucleotide.html
 * the DNAfull scoring matrix (uses IUPAC notation for ambiguous nucleotides)
 * gap opening penalty of 10 and gap extension penalty of 1
 * do not forget to count gaps on the ends using '-endweight'

Given: Two GenBank IDs.

Return: The maximum global alignment score between the DNA strings
associated with these IDs.

Sample Dataset
--------------

JX205496.1 JX469991.1

Sample Output
-------------

257

"""

import sys, re, os
from subprocess import check_output
from Bio import Entrez, SeqIO

"""
Solved using 'needle' from Debian 'emboss' package (6.6 in Jessie).
This Python script just reads the dataset file, downloads the sequences
from the Genbank (nuccore) database (needs a write access to the local
directory), feeds the sequence files to 'needle' with correct options
and parses the result. Finally, downloaded sequences are deleted. 
"""

# extract data from the input file
with open(sys.argv[1], 'r') as in_file:
  genbank_ids = in_file.read().split()

# specify your email address to stop warnings
Entrez.email = "kub4@mailinator.com"

# define the fasta filenames to use
filenames = ["a.fasta", "b.fasta"]

# create the database handle (see 03-FRMT for Entrez.efetch parameters)
handle = Entrez.efetch(db="nuccore", id=genbank_ids, rettype="fasta")

# parse using SeqIO and write both sequences to the local directory
for record, filename in zip(SeqIO.parse(handle, "fasta"), filenames):
  with open(filename, 'w') as out_file:
    SeqIO.write(record, out_file, "fasta")

# create a list for construction of the 'needle' command
command = []

# start with 'needle' followed by two fasta files
command.append("needle {} {}".format(filenames[0], filenames[1]))

# turn off the prompts and write the output to stdout
command.append("-auto -stdout")

# add the scoring matrix file (/usr/share/EMBOSS/data/EDNAFULL)
command.append("-datafile EDNAFULL")

# add gap opening penalty of 10 and gap extension penalty of 1
command.append("-gapopen 10 -gapextend 1")

# do not forget to count gaps on the ends
command.append("-endweight true -endopen 10 -endextend 1")

# construct the final command
command = ' '.join(command)

# run the command with arguments and return its output as a list
needleout = check_output(command, shell=True).decode().splitlines()

# create a regular expression for number recovery (int or float)
pattern = r'([-]?(([0-9]+(\.[0-9]+)?)|(\.[0-9]+)))'

# now scan the needleout and print the result
for line in needleout:
  if "# Score:" in line:
    print(str(float(re.search(pattern,line).group())))
    break

# remove the temporary sequence files
for filename in filenames:
  os.remove(filename)

# amen, kallisti
