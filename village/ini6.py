#!/usr/bin/env python3

"""
Dictionaries
============

Given: A string s of length at most 10000 letters.

Return: How many times any word occurred in the string.
Case sensitive! Lines in output can be in any order.

Sample Dataset
--------------

We tried list and we tried dicts also we tried Zen

Sample Output
-------------

and 1
We 1
tried 3
dicts 1
list 1
we 2
also 1
Zen 1

"""

import sys
from collections import Counter

with open(sys.argv[1], 'r') as in_file:   # to automagically close the file
  words = in_file.read().split()    # make a list of all words

c = Counter(words) # using high power magic, when available :)

for key, value in c.items():
  print (key, value)
