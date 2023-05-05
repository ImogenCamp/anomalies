#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 10:23:50 2023

@author: imogencamp
"""

"""This code has been written to produce solutions to the diophantine equation

        \sum_i n_i x_i^2 - \sum_j m_j y_j^2 = 0
    
subject to the constraint

        \sum_i n_i - \sum_j m_j = 0.

The initial solution used is

        [X] = [1:1:...:1].
        
We check that the code works using the explicit example

        x^2 + y^2 = 2z^2.
        
Note there is almost certainly a better way to remove the duplicate solutions...

"""

import numpy as np
from collections import Counter
import time

#we wish to time how long it takes the code to run
start_time = time.perf_counter()



"""useful mathematial functions below"""


#function to express x as a vector of integers whith a gcd of 1
def Divide_by_gcd(vector):
    vector = np.asarray(vector)
    sign = np.sign(vector)
    vector = np.abs(vector)
    gcd_value = np.gcd.reduce(vector.astype(int))
    if gcd_value != 0: 
        return (sign * vector / gcd_value).astype(int)
    else:
        return vector

#function to remove duplicate vectors in array of vectors
def Remove_duplicates(arrays):
    unique = []
    for vector in arrays:
        unique_vector = Divide_by_gcd(vector)
        if not any([np.array_equal(uv, unique_vector) for uv in unique]):
            unique.append(unique_vector)
    return np.array(unique)


#function to remove vectors which are integer multiples of one another in an array of vectors
def Remove_integer_multiples(arrays):
    unique_arrays = []
    arrays_set = set()
    for i, array1 in enumerate(arrays):
        abs_array1 = tuple(np.abs(Divide_by_gcd(array1)))
        if abs_array1 not in arrays_set:
            unique_arrays.append(array1)
            arrays_set.add(abs_array1)
    return np.array(unique_arrays)

    

"""useful functions for generating solutions below"""

#useful function for our solution parametrisation
def Psi_vector(coeffs, lambda_vector):
    squared_vector = np.power(lambda_vector, 2)
    return np.dot(coeffs, squared_vector)

#useful function for our solution parametrisation
def Phi_vector(coeffs, lambda_vector):
    return np.dot(coeffs, lambda_vector)


"""now genenerate the solution vector"""

#our parametrisation depends on psi, phi and lambda
#form derived in notes
def Solution_vector(psi_vector, phi_vector, lambda_vector):
    multiplied_vectors = 2 * phi_vector * lambda_vector
    return multiplied_vectors - psi_vector


"""now check that the diophantine equation is satisfied"""
#coefficients is a vector of the coefficients in the diophantine equation
#sol is a solution vector previously abtained
def Diophantine_checker(coeffs, sol):
    squared_sol = np.power(sol, 2)
    diophantine_sum = np.dot(coeffs, squared_sol)
    if diophantine_sum == 0:
        #print("True")
        return True
    else:
        #print("False")
        return False
    
def Mutual_diophantine_checker(coeffs, sol1, sol2):
    diophantine_sum = np.dot(coeffs, sol1 * sol2)
    if diophantine_sum == 0:
        #print("True")
        return True
    else:
        #print("False")
        return False
    
    
    
    
    
    
    
    

"""we now wish to check that this works for some examples"""

#define our inputs
coefficient_vector = np.array([1, 1, -2])
#check that input coefficients sum to zero (or our method of chords will not work)
if np.sum(coefficient_vector)!=0:
    raise ValueError("GRAVITATIONAL ANOMALY CANCELLATION NOT SATISFIED")

#define range of lambda values over which we want to scan
min_value = 0
max_value = 4000

#do this for a large number of possible lambda vectors
#always have final element vanishing so that plane does not contain x
#remember we are in projective space! should remove equivalent values of lambda (related by a scaling)
lambda_matrix = np.array([[a, b, 0] for a in range(min_value, max_value) for b in range(min_value, max_value)])
print("Number of lambda values is:", np.shape(lambda_matrix)[0])
lambda_matrix = Remove_integer_multiples(lambda_matrix)
lambda_matrix = np.array([array for array in lambda_matrix if not np.all(array == [0, 0, 0])])
print("Number of lambda values after duplicate removal is:",np.shape(lambda_matrix)[0])
checking_array=np.append(lambda_matrix, np.negative(lambda_matrix), axis = 0)


#check for duplicate lambda values

list_of_tuples = [tuple(x) for x in checking_array] #list of solutions, checking array added to make sure that solutions related by negative multiples have also been removed
counts = Counter(list_of_tuples) # count occurrences of each tuple
duplicates = [tuple for tuple, count in counts.items() if count > 1] # get tuples with count greater than 1

if duplicates:
    print("There ARE duplicate lambda values:", duplicates)
else:
    print("There AREN'T duplicate lambda values.")
    
counts = {}

for t in list_of_tuples:
    if t in counts:
        counts[t] += 1
    else:
        counts[t] = 1
        
for t, count in counts.items():
    if count > 1:
        print("The tuple", t, "appears", count, "times.")



#initialize a numpy array of booleans
diophantine_boolean_array = []
#initialize an array for solution vectors
#convert this to numpy array *after* looping (apparently this is faster)
diophantine_solution_array = []

  
for i in range(0, int(lambda_matrix.shape[0])):
    lambda_vector = lambda_matrix[i]
    
    #define our useful variables
    psi_vector = Psi_vector(coefficient_vector, lambda_vector)
    phi_vector = Phi_vector(coefficient_vector, lambda_vector)

    #output solution
    #remember we are in projective space, so all solutions related by scaling are equivalent
    #choose to output solution with integers with gcd of 1
    solution_vector = Divide_by_gcd(Solution_vector(psi_vector, phi_vector, lambda_vector))
    diophantine_solution_array.append(solution_vector)
    
    #print solution and check it satisfies the diophantine equation
    #diophantine_boolean_array = np.append(diophantine_boolean_array, Diophantine_checker(coefficient_vector, solution_vector))
    diophantine_boolean_array.append(Diophantine_checker(coefficient_vector, solution_vector))
 




print("Number of solutions is:", np.shape(lambda_matrix)[0])
    
#remove duplicate solutions
diophantine_solution_array = Remove_integer_multiples(diophantine_solution_array)
diophantine_solution_array = [array for array in diophantine_solution_array if not np.all(array == [0, 0, 0])]


print("Number of solutions after duplicate removal is:",np.shape(lambda_matrix)[0])

#check for any duplicate solutions (i.e. did this work?)

checking_array_2 = np.append(diophantine_solution_array, np.negative(diophantine_solution_array), axis = 0)
list_of_tuples = [tuple(x) for x in checking_array_2] #list of solutions
counts = Counter(list_of_tuples) # count occurrences of each tuple
duplicates = [tuple for tuple, count in counts.items() if count > 1] # get tuples with count greater than 1

if duplicates:
    print("There ARE duplicate solutions:", duplicates)
else:
    print("There AREN'T duplicate solutions.")
    
counts = {}

for t in list_of_tuples:
    if t in counts:
        counts[t] += 1
    else:
        counts[t] = 1
        
for t, count in counts.items():
    if count > 1:
        print("The tuple", t, "appears", count, "times.")
        
#there should be no duplicate solutions

 
 


#now convert our arrays to numpy arrays (i think this will make them easier to deal with, although i am not sure)
diophantine_solution_array = np.asarray(diophantine_solution_array)
diophantine_boolean_array = np.asarray(diophantine_boolean_array)


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

"""
output

    Number of lambda values is: 16000000
    Number of lambda values after duplicate removal is: 9724005
    There AREN'T duplicate lambda values.
    Number of solutions is: 9724005
    Number of solutions after duplicate removal is: 9724005
    There AREN'T duplicate solutions.
    YAY IT WORKS
    Elapsed time: 470.3239 seconds
"""












