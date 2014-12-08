#!/usr/bin/env python3

"""
Reverse Complement Problem
==========================


Sample Dataset
--------------

AAAACCCGGT

Sample Output
-------------

ACCGGGTTTT

"""
import sys

with open(sys.argv[1], 'r') as input_file:
  input_data = input_file.read()

