#!/usr/bin/env python3

"""
Rabbits and Recurrence Relations
================================

A recurrence relation is a way of defining the terms of a sequence with respect
to the values of previous terms. In the case of Fibonacci's rabbits, any given
month will contain the rabbits that were alive the previous month, plus any
new offspring. A key observation is that the number of offspring in any month
is equal to the number of rabbits that were alive two months prior. As a
result, if Fn represents the number of rabbit pairs alive after the n-th month,
then we obtain the Fibonacci sequence having terms Fn that are defined by the
recurrence relation Fn=F(n−1)+F(n−2) (with F1=F2=1 to initiate the sequence).

When finding the n-th term of a sequence defined by a recurrence relation, we
can simply use the recurrence relation to generate terms for progressively
larger values of n. This problem introduces us to the computational technique
of dynamic programming, which successively builds up solutions by using the
answers to smaller cases.

Given: Positive integers n≤40 and k≤5.

Return: The total number of rabbit pairs that will be present after n months
if we begin with 1 pair and in each generation, every pair of reproduction-age
rabbits produces a litter of k rabbit pairs (instead of only 1 pair).

Sample Dataset
--------------

5 3

Sample Output
-------------

19

"""
import sys

with open(sys.argv[1], 'r') as input_file:   # to automagically close the file
  input_data = input_file.read()             # when leaving the nested block

input_data = input_data.split()
n = int(input_data[0])
k = int(input_data[1])

assert (n > 0), "n should be positive!"
assert (k > 0), "k should be positive!"

rabbits = [1, 1]                            # F1=F2=1 to initiate the sequence

for i in range(2, n):                       # this is okay even for n=1 :)
  new_rabbits = rabbits[i-1] + k * rabbits[i-2]
  rabbits.append(new_rabbits)

print(rabbits[n-1])
