#!/usr/bin/env python3
"""
Pairwise Local Alignment
========================

A local alignment of two strings s and t is an alignment of substrings r and
u of s and t, respectively. Finding an optimal local alignment involves two
maximizations: each pair of substrings has a maximum score (for their global
alignment), and we want the pair of substrings with the largest possible
maximum score. 

Use:
 * http://www.ebi.ac.uk/Tools/psa/emboss_water/
 * the BLOSUM62 scoring matrix
 * gap opening penalty of 10
 * gap extension penalty of 1

Given: Two UniProt ID's corresponding to two protein strings s and t.

Return: The maximum score of any local alignment of s and t. 

Sample Dataset
--------------

B3ET80 Q78PG9

Sample Output
-------------

35

"""

import sys, re, os
from Bio import ExPASy
from Bio import SwissProt
from subprocess import check_output

"""
Solved using 'water' from Debian 'emboss' package (6.6 in Jessie).
This Python script just reads the dataset file, downloads the sequences
from the ExPASy database (needs a write access to the local directory),
feeds the sequence files to 'water' with correct arguments and parses
the result. Finally, downloaded sequences are deleted.

NOTE: the answer is converted to int using the round function, as it seems
that the grader does not accept any floats. With the penalties given, this
is only matter of formatting, but with different penalties, the result
might be become inaccurate !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
"""

# extract data from the input file
with open(sys.argv[1], "r") as in_file:
  uniprot_ids = in_file.read().split()

# creating one handle for more proteins does not seem to be supported
# (the comment in the rosalind's dbpr programming shortcut is wrong,
# if I understand the surprisingly simple biopython's source correctly)
for pid in uniprot_ids:
  record = SwissProt.read(ExPASy.get_sprot_raw(pid))
  with open(pid, "w") as out_file:      # using pid as a filename
    out_file.write(">" + pid + "\n")    # fasta header
    out_file.write(record.sequence)     # sequence

# create a list for construction of the 'water' command
command = []

# start with 'water' followed by two fasta files
command.append("water {} {}".format(uniprot_ids[0], uniprot_ids[1]))

# turn off the prompts and write the output to stdout
command.append("-auto -stdout")

# add the scoring matrix file (/usr/share/EMBOSS/data/EBLOSUM62)
command.append("-datafile EBLOSUM62 ")

# add gap opening penalty of 10 and gap extension penalty of 1
command.append("-gapopen 10 -gapextend 1")

# construct the final command
command = ' '.join(command)

# run the command with arguments and return its output as a list
waterout = check_output(command, shell=True).decode().splitlines()

# create a regular expression for number recovery (int or float)
pattern = r'([-]?(([0-9]+(\.[0-9]+)?)|(\.[0-9]+)))'

# now scan the waterout and print the result
for line in waterout:
  if "# Score:" in line:
    answer = float(re.search(pattern,line).group())
    answer = round(answer) #for some reason, the grader accepts only ints!!!
    print(answer)
    break

# remove the temporary sequence files
for filename in uniprot_ids:
  os.remove(filename)
