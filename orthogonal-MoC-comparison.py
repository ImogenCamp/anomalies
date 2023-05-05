#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 11:07:42 2023

@author: imogencamp
"""
"""The goal of this code is to compare the solution sets produced by the Method of Chords and the generalised
orthogonal matrix method. We see that the intersection between two sets of size O(10^5) is approximately 10%.
Not too much can be concluded from this, as the solution sets are parameterised by different parameters. There
are, after all, infinitely many rational numbers.
"""



import numpy as np
import time
from fractions import Fraction

#we wish to time how long it takes the code to run
start_time = time.perf_counter()


"""define functions"""


#function to remove duplicate vectors in array of vectors
def Remove_duplicates(arrays):
    unique = []
    for vector in arrays:
        unique_vector = Divide_by_gcd(vector)
        if not any([np.array_equal(uv, unique_vector) for uv in unique]):
            unique.append(unique_vector)
    return np.array(unique)


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
    #if diophantine_sum == 0:
    if np.isclose(diophantine_sum, 0, rtol=1e-6, atol=1e-6):
        return True
    else:
        return False
    
#define a function to check orthogonality
def Orthogonal_checker(matrix, B):
    #tolerence chosen to be 0.005 because of errors
    if np.allclose(matrix.T @ B @ matrix, B, rtol=0.005):
        return True
    else:
        return False
   
#function to convert floats to fractions
def to_fraction(x):
    return Fraction.from_float(x).limit_denominator()


#function converts each float in an array to a fraction
vfunc = np.vectorize(to_fraction)
    

"""start by producing an array of solutions from the method of chords"""

#define range to scan over
min_value = 1
max_value = 100

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







"""now produce solutions using the orthogonal matrices"""

#define quadratic form matrix, its inverse, and the identity
B = np.diag([1, 1, -1, -3]) 
B_inv = np.linalg.inv(B)
I = np.identity(4)


#define an array of matrices that I will scan over
matrices = []

#define range of values to scan over
min_value = 1
max_value = 10

for a in range(min_value, max_value):
    for b in range(min_value, max_value):
        for c in range(min_value, max_value):
            for d in range(min_value, max_value):
                for e in range(min_value, max_value):
                    for f in range(min_value, max_value):
                        #define matrix according to constraints in Workings 25
                        A = np.array([[0, a, b, c], [-a, 0, d, e], [b, d, 0, f], [c/3, e/3, -f/3, 0]])
                        
                        #compute the eigenvalues of A
                        eigenvalues = np.linalg.eigvals(A)
                        
                        #check if any of the eigenvalues are -1
                        is_negative_one = np.any(np.isclose(eigenvalues, -1))
                        
                        #we do not want to consider matrices with -1 as as an eigenvalue, since eigenvalue of (A + I) is 0, cannot take inverse
                        if not is_negative_one:
                            matrices.append(A)
                            

"""now use Cayley parameterisation to find orthogonal matrices"""

#array for orthogonal matrices
orthogonal_matrices = []
#array to check that they are orthogonal
orthogonal_boolean_array = []

for i in range(len(matrices)):
    #for each matrix we added to matrices
    A = matrices[i]
    
    #define the orthogonal matrix (see Cayley parameterisation) and check if orthogonal
    O = np.linalg.inv(I + A) @ (I - A)
    orthogonal_bool = Orthogonal_checker(O, B)
    
    #add O to solution array and orthogonal_bool to boolean array
    orthogonal_matrices.append(O)
    orthogonal_boolean_array.append(orthogonal_bool)


if np.all(orthogonal_boolean_array)  == True:
    print("YAY IT WORKS")
else:
    print("O no")
    
    
#define initial solution
vector = np.array([1, 0, 1, 0])
#array of new solutions is:
vector_array = np.array([np.dot(matrix, vector) for matrix in orthogonal_matrices])

checker_array = np.array([Diophantine_checker(np.array([1, 1, -1, -3]), sol) for sol in vector_array])


if np.all(checker_array):
    print("YAY IT WORKS")
else:
    print("O no")


orthogonal_solution_array = []


#convert to integers for comparison
for i in range(0, vector_array.shape[0]):
    sol_vec = vector_array[i]
    #remove zero components (denominator will be undefined otherwise)
    sol_vec_nonzero = [flt for flt in sol_vec if not np.isclose(flt, 0, rtol=1e-6, atol=1e-6)]
    #convert to an array of fractions
    as_fractions = vfunc(sol_vec_nonzero)
    # convert fractions to denominators
    denominators = np.array([f.denominator for f in as_fractions])
    # find LCM of denominators
    lcm = np.lcm.reduce(denominators)
    lcm = float(lcm)
    #integer solution array!
    sol_vec = np.array([np.round(lcm * x).astype(int) for x in sol_vec])
    #append the integer solution array to out final array
    orthogonal_solution_array.append(sol_vec)


#print(orthogonal_solution_array)
L = len(orthogonal_solution_array)
print(L)


checker_array_2 = np.array([Diophantine_checker(np.array([1, 1, -1, -3]), orthogonal_solution_array[i]) for i in range(0, L)])


if np.all(checker_array_2):
    print("YAY IT WORKS")
else:
    print("O no")





"""now we need to compare orthogonal_solution_array and diophantine_solution_array for their overlap"""


#orthogonal_solution_array = np.array([orthogonal_solution_array])
orthogonal_solution_array = tuple(map(tuple, orthogonal_solution_array))
#diophantine_solution_array = np.array([diophantine_solution_array])
diophantine_solution_array = tuple(map(tuple, diophantine_solution_array))

#define sets to remove duplicate elements and so we can define 
orthogonal_set = set(tup for tup in orthogonal_solution_array)
diophantine_set = set(tup for tup in diophantine_solution_array)

print("number of solutions from matrices is:", len(orthogonal_set))
print("number of solutions from MoC is:", len(diophantine_set))

intersection = orthogonal_set.intersection(diophantine_set)
print("size of intersection is:", len(intersection))


#print(intersection)
#intersection is quite small...

#print(orthogonal_set - intersection)






"""runtime"""

#output time taken for code to run
end_time = time.perf_counter()
elapsed_time = end_time - start_time

print(f"Elapsed time: {elapsed_time:.4f} seconds")

"""Output:
    
    970299
    YAY IT WORKS
    YAY IT WORKS
    YAY IT WORKS
    529766
    YAY IT WORKS
    number of solutions from matrices is: 105492
    number of solutions from MoC is: 811213
    size of intersection is: 7185
    Elapsed time: 169.1979 seconds
    
"""








