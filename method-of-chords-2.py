#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 18:07:10 2023

@author: imogencamp
"""


"""This is an example of using the method of chords in a Witt decomposed system

        x**2 + y**2 - w**2 - 3z**2 = 0.
        
Why? So we can more easily compare our solutions from the orthogonal matrices without having to 
transform back to the original basis.

Matrix defining quadratic form:
    
        B = np.diag(1, 1, -1, -3).

Initial solution:
    
        [X] = [1 : 0 : 1 : 0].
        
Parameters:
    
        [A] = [0 : a1 : a2 : a3].
        
Solution:
    
        [Y] = [3*a3**2 + a2**2 - a1 : 2*a1*a2 : 3*a3**2 + 3*a2**2 - a1 : 2*a3*a2].
        
"""

import numpy as np
import time

#we wish to time how long it takes the code to run
start_time = time.perf_counter()


#define function to express x as a vector of integers whith a gcd of 1
def Divide_by_gcd(vector):
    vector = np.asarray(vector)
    sign = np.sign(vector)
    vector = np.abs(vector)
    gcd_value = np.gcd.reduce(vector.astype(int))
    if gcd_value != 0: 
        return (sign * vector / gcd_value).astype(int)
    else:
        return vector

#define function to compute solutions
def Solution_method_of_chords(arr):
    a1 = arr[0]
    a2 = arr[1]
    a3 = arr[2]
    solution_array = np.array([a1**2 - a2**2 - 3*a3**2, 2*a1*a2, a1**2 + a2**2 - 3*a3**2, 2*a3*a2])
    return Divide_by_gcd(solution_array)


#function to check if diiophantine equation is satisfied
def Diophantine_checker(B, sol):
    squared_sol = np.power(sol, 2)
    diophantine_sum = np.dot(B, squared_sol)
    if diophantine_sum == 0:
        return True
    else:
        return False
    

"""Now perform checks."""

#define range to scan over
min_value = 1
max_value = 10

B = np.array([1, 1, -1, -3])
lambda_matrix = np.array([[a, b, c] for a in range(min_value, max_value) for b in range(min_value, max_value) for c in range(min_value, max_value)])



#initialize a numpy array of booleans
diophantine_boolean_array = []
#initialize an array for solution vectors
diophantine_solution_array = []

for i in range(0, int(lambda_matrix.shape[0])):
    #output solution
    solution_vector = Solution_method_of_chords(lambda_matrix[i])
    
    #append to arrays
    diophantine_solution_array.append(solution_vector)
    diophantine_boolean_array.append(Diophantine_checker(B, solution_vector))
    
print(len(diophantine_solution_array))
    
 
#check that all our solutions are indeed solutions...
if np.all(diophantine_boolean_array)  == True:
    print("YAY IT WORKS")
    #print(diophantine_solution_array)
else:
    print("O no")



#output time taken for code to run
end_time = time.perf_counter()
elapsed_time = end_time - start_time

print(f"Elapsed time: {elapsed_time:.4f} seconds")
















