#!/usr/bin/env python3

"""
Fibonacci Numbers
=================

The Fibonacci numbers 0,1,1,2,3,5,8,13,21,34,…

Given: A positive integer n≤25.

Return: The value of Fn.

Sample Dataset
--------------

6

Sample Output
-------------

8

"""
import sys

# define Fibonacci function
def fib(n):
  if n==0:
    return 0
  f=[0, 1]
  for i in range(2, n+1):
    f.append(f[i-1]+f[i-2])
  return f[n]

"""
#the above function is much faster (polynomial)
#than the naive recursive approach (exponential): 
"""
# define Fibonacci function, the naive way
def fib_naive(n):
  if n==0:
    return 0
  if n==1:
    return 1
  return fib_naive(n-1) + fib_naive(n-2)

# use the command line argument as n (if int)
# or use it to open a file, then sanity check
try:
  n = int(sys.argv[1])
except:
  with open(sys.argv[1], 'r') as in_file:
    n = int(in_file.read().strip())

assert n>=0, "Need positive n!"

# call Fibonacci function
print(fib(n))
#print(fib_naive(n))
