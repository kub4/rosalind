#!/usr/bin/env python3
"""
Count programming languages used by Rosalind's TOP100 users. Outputs both
plain stats and weighted stats (100 points for the first user, one less
for every next place.
"""

from urllib.request import urlopen
from collections import defaultdict
import datetime

top100url = "http://rosalind.info/statistics/top/"
top100users = []

with urlopen(top100url) as top100:
  for line in top100.readlines():
    line = line.decode().strip()
    if 'href="/users' in line:
      top100users.append(line)

for i, user in enumerate(top100users):
   start = user.find('href="/users')
   end = user.find('">')
   top100users[i] = user[start+13:end-1]

today = datetime.date.today()
print("Rosalind.info TOP100 users", today)

for i, user in enumerate(top100users):
  url = "http://rosalind.info/users/"+user+"/"
  print(url)
  with urlopen(url) as homepage:
    lang = 0
    language = None
    for line in homepage.readlines():
      line = line.decode().strip()
      if lang == 1:
        language = line
        break
      if "<dt>Language</dt>" in line:
        lang = 1
  if language:
    start = language.find('<dd>')
    stop = language.find('</dd>')
    language = language[start+4:stop]
  top100users[i] = [top100users[i], language]

languages = defaultdict(int)
languages_weighted = defaultdict(int)

for i, userdata in enumerate(top100users):
  lang = userdata[1]
  languages[lang] += 1
  languages_weighted[lang] += (100-i)

print("\nLANGUAGES:")
for key, value in sorted(languages.items(),key=lambda kv: kv[1], reverse=True):
  print(key, value)

print("\nLANGUAGES WEIGHTED:")
for key, value in sorted(languages_weighted.items(),key=lambda kv: kv[1], reverse=True):
  print(key, value)
