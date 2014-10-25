#!/usr/bin/env python3

"""
Mortal Fibonacci Rabbits
========================

Recall the definition of the Fibonacci numbers from “Rabbits and Recurrence
Relations”, which followed the recurrence relation Fn=F(n−1)+F(n−2) and assumed
that each pair of rabbits reaches maturity in one month and produces a single
pair of offspring (one male, one female) each subsequent month.

Our aim is to somehow modify this recurrence relation to achieve a dynamic
programming solution in the case that all rabbits die out after a fixed number
of months. See Figure 4 for a depiction of a rabbit tree in which rabbits live
for three months (meaning that they reproduce only twice before dying).

Given: Positive integers n≤100 and m≤20.

Return: The total number of pairs of rabbits that will remain after the n-th
month if all rabbits live for m months.

Sample Dataset
--------------

6 3

Sample Output
-------------

4

"""
import sys

with open(sys.argv[1], 'r') as input_file:   # to automagically close the file
  input_data = input_file.read()             # when leaving the nested block

input_data = input_data.split()
n = int(input_data[0])
m = int(input_data[1])

assert (n > 0), "n should be positive!"
assert (m > 0), "m should be positive!"

##############################################################################
# well, this is probably very naive and inefficient, but easy to understand, #
# so we are going to store the complete age structure for each generation:   #
##############################################################################

gens = []                         # creating a list of all generations (empty)

first_gen = [1]+[0]*(m-1)         # initializing the first generation
gens.append(first_gen)            # append to the list of generations

for g in range(1, n):             # for each next generation...
  # sum all adults of the previous gen, use as newborns, and copy the rest
  # from the previous gen, shifted by one age group, discarding the oldest
  new_gen = [sum(gens[g-1][1:])] + gens[g-1][:-1]
  gens.append(new_gen)            # append to the list of generations

print(sum(gens[n-1]))
