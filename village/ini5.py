#!/usr/bin/env python3

"""
Working with Files
==================

Given: A file containing at most 1000 lines.

Return: A file containing all the even-numbered lines from the original file.
Assume 1-based numbering of lines.

Sample Dataset
--------------

Bravely bold Sir Robin rode forth from Camelot
Yes, brave Sir Robin turned about
He was not afraid to die, O brave Sir Robin
And gallantly he chickened out
He was not at all afraid to be killed in nasty ways
Bravely talking to his feet
Brave, brave, brave, brave Sir Robin
He beat a very brave retreat

Sample Output
-------------

Yes, brave Sir Robin turned about
And gallantly he chickened out
Bravely talking to his feet
He beat a very brave retreat

"""

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
