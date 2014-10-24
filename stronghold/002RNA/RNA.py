#!/usr/bin/env python3

"""
========================


"""

import sys

with open(sys.argv[1], 'r') as input_file:   # to automagically close the file
  data = input_file.read()                   # when leaving the nested block

