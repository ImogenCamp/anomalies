#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 12:29:27 2023

@author: imogencamp
"""

"""This is an example of using the generalised orthogonal matrix method in a Witt decomposed system

        x**2 + y**2 - w**2 - 3z**2 = 0.

Matrix defining quadratic form:
    
        B = np.diag(1, 1, -1, -3).

Initial solution:
    
        [X] = [1 : 0 : 1 : 0].
        
Generalised skew-symmetric matrix:
    
        A = {{0, a, b, c}, {-a, 0, d, e}, {b, d, 0, f}, {c/3, e/3, -f/3, 0}}
        
It may very well be quicker to find the generalised orthogonal matrices and multiply the initial solution
in Mathematica. Then we would just have to input the solution and scan over the variables. However, this 
is a bit of a pain to do (differing syntax, long solution...), and so I have not.
        
"""

import numpy as np
import time

#we wish to time how long it takes the code to run
start_time = time.perf_counter()

#define a function to check orthogonality
def Orthogonal_checker(matrix, B):
    #tolerence chosen to be 0.005 because of errors
    if np.allclose(matrix.T @ B @ matrix, B, rtol=0.005):
        return True
    else:
        return False


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
                            
print(len(matrices))

"""now use Cayley parameterisation to find orthogonal matrices"""

#array for orthogonal matrices
orthogonal_matrices = []
#array to check that they are orthogonal
boolean_array = []

for i in range(len(matrices)):
    #for each matrix we added to matrices
    A = matrices[i]
    
    #define the orthogonal matrix (see Cayley parameterisation) and check if orthogonal
    O = np.linalg.inv(I + A) @ (I - A)
    orthogonal_bool = Orthogonal_checker(O, B)
    
    #add O to solution array and orthogonal_bool to boolean array
    orthogonal_matrices.append(O)
    boolean_array.append(orthogonal_bool)


if np.all(boolean_array)  == True:
    print("YAY IT WORKS")
else:
    print("O no")



#output time taken for code to run
end_time = time.perf_counter()
elapsed_time = end_time - start_time

print(f"Elapsed time: {elapsed_time:.4f} seconds")


"""Output:
    
    529766
    YAY IT WORKS
    Elapsed time: 44.7287 seconds
    
"""


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    