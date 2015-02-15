#!/usr/bin/env python3

"""
Independent Segregation of Chromosomes
======================================

Consider a collection of coin flips. One of the most natural questions we can
ask is if we flip a coin 92 times, what is the probability of obtaining 51
"heads", vs. 27 "heads", vs. 92 "heads"?

Each coin flip can be modeled by a uniform random variable in which each of the
two outcomes ("heads" and "tails") has probability equal to 1/2. We may assume
that these random variables are independent; in layman's terms, the outcomes
of the two coin flips do not influence each other.

A binomial random variable X takes a value of k if n consecutive "coin flips"
result in k total "heads" and n−k total "tails." We write that X∈Bin(n,1/2).

Given: A positive integer n≤50.

Return: An array A of length 2n in which A[k] represents the common logarithm
of the probability that two diploid siblings share at least k of their 2n
chromosomes (we do not consider recombination for now).

Sample Dataset
--------------

5

Sample Output
-------------

0.000 -0.004 -0.024 -0.082 -0.206 -0.424 -0.765 -1.262 -1.969 -3.010

"""

import sys
from math import log
from math import factorial as fac

# define a function for binomial probability mass function
# (the probability of getting _exactly_ k successes in n trials)
# http://en.wikipedia.org/wiki/Binomial_distribution
def binomial_prob_mass(k,n,p):
  return (fac(n)/(fac(k)*fac(n-k)))*(p**k)*((1-p)**(n-k))

# open a file and get n (in the genetics sense)
with open(sys.argv[1], 'r') as in_file:
  gen_n = int(in_file.read().strip())

# what we get from input is n in the genetics sense
# but we need n in the combinatorics sense which
# is equal to the total number of chromosomes (2n)
n = 2 * gen_n

# a list to keep the binomial probabilities (non-cumulative)
probs = []

# it is not very clear from the problem description, but after
# some trial and error I know that they are not interested in
# k=0, so we will use k from 1 to n...
for k in range(1, n+1):
  probs.append(binomial_prob_mass(k, n, 0.5))

# now compute the cumulative probabilities by summing
# the non-cumulative probabilities...
cumprobs = []

for i in range(0, len(probs)):
  cumprobs.append(sum(probs[i:])) 

# print the results in log form, 3 decimal places
# there are some differences from the sample output due to
# precision issues (may be even on the rosalind side),
# but rosalind allows +/- 0.001 inaccuracy, so no stress... 
print(" ".join(["{:.3f}".format(log(x,10)) for x in cumprobs]))
