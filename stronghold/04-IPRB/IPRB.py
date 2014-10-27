#!/usr/bin/env python3

"""
Mendel's First Law
==================

Given: Three positive integers k, m, and n, representing a population
containing k+m+n organisms: k individuals are homozygous dominant for a factor,
m are heterozygous, and n are homozygous recessive.

Return: The probability that two randomly selected mating organisms will
produce an individual possessing a dominant allele (and thus displaying the
dominant phenotype). Assume that any two organisms can mate.

Sample Dataset
--------------

2 2 2

Sample Output
-------------

0.78333

"""
import sys

with open(sys.argv[1], 'r') as in_file:   # to automagically close the file
  data = list(map(int,in_file.readline().split()))

k = data[0]   # dominant homozygotes (AA)
m = data[1]   # heterozygotes (Aa)
n = data[2]   # recessive homozygotes (aa)
s = k+m+n     # sum of all individuals

assert (k > 0), "k should be positive!"
assert (m > 0), "m should be positive!"
assert (n > 0), "n should be positive!"
assert (s > 1), "need a pair to have sex!"

'''
a rewrite using a closed function, computing probability of recessive offspring
(and then substracting  the result from one):

let's call the first parent "father", the second "mother" (but of course, we
are dealing with hermafrodites)

1. if "father" is AA, there is no chance for a recessive child -> ignoring

2. if "father" is aa (n/s), all depends on mother:

 2a) "father" aa, "mother" AA -> no recessive children, ignoring

 2b) "father" aa, "mother" aa -> 100% of recessive children
     (n/s)*((n-1)/(s-1))
     -------------------
 
 2c) "father" aa, "mother" Aa -> 50% of recessive children
     (n/s)*(m/(s-1))*0.5
     -------------------

3. if "father" is Aa (m/s), mother is important:

 3a) "father" Aa, "mother" AA -> no recessive children, ignoring

 3b) "father" Aa, "mother" aa -> 50% of recessive children
     (m/s)*(n/(s-1))*0.5
     -------------------
 
 3c) "father" Aa, "mother" Aa -> 25% of recessive children
     (m/s)*((m-1)/(s-1))*0.25
     ------------------------

now we need to add all the underlined expressions

(n/s)*((n-1)/(s-1)) + (n/s)*(m/(s-1))*0.5 +
+ (m/s)*(n/(s-1))*0.5 + (m/s)*((m-1)/(s-1))*0.25

then simplify and substract from one...

'''

r = (n*(n-1) + n*m + 0.25*m*(m-1)) / (s*(s-1))

dom_child_probability = 1 - r

print(dom_child_probability)
