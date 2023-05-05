#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  4 15:22:17 2023

@author: imogencamp
"""
"""This code performs a brute force search for solutions to the Diophantine equation

        x^2 + y^2 = 2z^2.
        
The goal is to compare the efficiency of this scan with the parameterisation of the Method of Chords.
"""


import time

#we wish to time how long it takes the code to run
begin_time = time.perf_counter()

def diophantine(a, b, c, max_value, max_time):
    solutions = [] #i think numpy arrays are not so good for for loops, so define a normal array
    start_time = time.time()
    for x in range(1, max_value):
        for y in range(1, max_value):
            for z in range(1, max_value):
                if a * x * x + b * y * y == c * z * z:
                    solutions.append([x, y, z]) #append the solutions to the array
                    if time.time() - start_time >= max_time:
                        return solutions #code stops running after a given time 
    return solutions


a = 1
b = 1
c = 2
max_value = 10000
max_time = 480 #seconds


result = diophantine(a, b, c, max_value, max_time)
print(len(result))
#print(result)

#output time taken for code to run
end_time = time.perf_counter()
elapsed_time = end_time - begin_time

print(f"Elapsed time: {elapsed_time:.4f} seconds")

"""Output:
    
    191
    Elapsed time: 481.8865 seconds
"""