#!/usr/bin/env python3

"""
Insertion Sort
==============

Although it is one of the elementary sorting algorithms with O(n^2) worst-case
time, insertion sort is the algorithm of choice either when the data is nearly
sorted (because it is adaptive) or when the problem size is small (because it
has low overhead).

For these reasons, and because it is also stable, insertion sort is often used
as the recursive base case (when the problem size is small) for higher overhead
divide-and-conquer sorting algorithms, such as “Merge Sort” or “Quick Sort”.

Given: A positive integer n≤103 and an array A[1..n] of integers.

Return: The number of swaps performed by insertion sort algorithm on A[1..n].

Sample Dataset
--------------

6
6 10 4 5 1 2

Sample Output
-------------

12

"""
import sys

# define the insertion sort function
def insertion_sort(a):
  """
  Sorts the array (list) 'a' in place.
  Returns number of swaps performed.
  http://en.wikipedia.org/wiki/Insertion_sort
  """
  swaps = 0
  for i in range(1,len(a)):
    k = i
    while k>0 and a[k-1]>a[k]:
      a[k-1],a[k] = a[k],a[k-1]
      k = k-1
      swaps += 1
  return swaps

# open the input file and get the array
with open(sys.argv[1], 'r') as in_file:
  array = [int(x) for x in in_file.readlines()[1].split()]

# sort the array
swaps = insertion_sort(array)
#print(array)

# print the number of swaps performed
print(swaps)
