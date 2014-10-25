#!/usr/bin/env python3

import sys
from itertools import islice

with open(sys.argv[1], 'r') as in_file:      # to automagically close the file
  for line in islice(in_file, 1, None, 2):   # when leaving the nested block
    print(line, end="")

"""
islice(in_file, 1, None, 2) skips the first line (start=1), then iterates over
all lines (stop=None), returning every second line (step=2). This'll work with
whatever filesize you throw at it; it won't need any more memory than required
by the internal iterator buffer.
"""
