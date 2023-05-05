#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 10:23:50 2023

@author: imogencamp
"""

"""This piece of code is similar to that in method-of-chords.py, except we are now taking pairs of solutions
to the n=1 case and checking whether they satisfy an n=2 system. Note that here it is especially important 
that we remove any duplicates, as otherwise we will conclude linearly independent n=2 solutions exist when
they may very well not.
"""



import numpy as np
import time

#we wish to time how long it takes the code to run
start_time = time.perf_counter()



"""useful mathematial functions below
these are not functions i am 100% happy with, need to check"""


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

#function to remove vectors from an array that are related by a scaling
#actually this does the same thing as the above function, but my code works and I don't want to mess with it
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
#check that input sum to zero (or our method of chords will not work)
if np.sum(coefficient_vector)!=0:
    raise ValueError("GRAVITATIONAL ANOMALY CANCELLATION NOT SATISFIED")

#define range of lambda values over which we want to scan
min_value = 0
max_value = 300

#do this for a large number of possible lambda vectors
#always have final element vanishing so that plane does not contain x
"""remember we are in projective space! should remove equivalent values of lambda (related by a scaling)"""
lambda_matrix = np.array([[a, b, 0] for a in range(min_value, max_value) for b in range(min_value, max_value)])
lambda_matrix = Remove_integer_multiples(lambda_matrix)
lambda_matrix = np.array([array for array in lambda_matrix if not np.all(array == [0, 0, 0])])
print("Number of lambda values after duplicate removal is:",np.shape(lambda_matrix)[0])




#initialize a numpy array of booleans
#diophantine_boolean_array = np.array([], dtype=bool)
diophantine_boolean_array = []
#initialize an array for solution vectors
#convert this to numpy array *after* looping (apparently this is faster)
diophantine_solution_array = []

  
for i in range(0, int(lambda_matrix.shape[0])):
    lambda_vector = lambda_matrix[i]
    
    #define out useful variables
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
 
    
#remove duplicate solutions
diophantine_solution_array = Remove_integer_multiples(diophantine_solution_array)
diophantine_solution_array = [array for array in diophantine_solution_array if not np.all(array == [0, 0, 0])]


print("Number of solutions after duplicate removal is:",np.shape(lambda_matrix)[0])




#now convert our arrays to numpy arrays (i think this will make them easier to deal with, although i am not sure)
diophantine_solution_array = np.asarray(diophantine_solution_array)
diophantine_boolean_array = np.asarray(diophantine_boolean_array)


#check that all our solutions are indeed solutions...
if np.all(diophantine_boolean_array)  == True:
    print("YAY IT WORKS")
    #print(diophantine_solution_array)
else:
    print("O no")

#YAY IT WORKS



"""Now we are going to test if any pairs of these solutions 
can be taken as solutions for U(1)^2."""


#find number of elements that we need to iterate over
rows = diophantine_solution_array.shape[0]

#initialise a flag to indicate whether any PAIR of solutions satisfy the system of diophantine equations
#these pairs will be called "mutual solutions" due to my lacklustre vocabulary
has_pair = False
#initialise mutual solution array
mutual_solution_array = []

#we can put any pair of solutions into Mutual_Diophantine_Checker
for i in range(rows):
    sol1 = diophantine_solution_array[i]
    for j in range(i + 1, rows):
        sol2 = diophantine_solution_array[j]
        if Mutual_diophantine_checker(coefficient_vector, sol1, sol2) == True:
            has_pair = True
            list1 = [list(diophantine_solution_array[i]), list(diophantine_solution_array[j])]
            mutual_solution_array.append(list1)


if has_pair == False:
    print("no U(1)^2 solutions")
else:
    print("solutions exist fr U(1)^2")
    print(mutual_solution_array)


#output time taken for code to run
end_time = time.perf_counter()
elapsed_time = end_time - start_time

print(f"Elapsed time: {elapsed_time:.4f} seconds")

"""output:
    
    Number of lambda values after duplicate removal is: 54637
    Number of solutions after duplicate removal is: 54637
    YAY IT WORKS
    no U(1)^2 solutions
    Elapsed time: 976.8157 seconds
"""










