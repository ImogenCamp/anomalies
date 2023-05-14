#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  8 00:05:22 2023

@author: imogencamp
"""


"""This piece of code is written to implement the algorithm I developed to explicitly perform a 
Witt decomposition. Note that it is necessary to import functions from the file witt_index.py."""

import numpy as np
import time
from sympy import Matrix
from witt_index import *
from fractions import Fraction
from sympy import factorint


#we wish to time how long it takes the code to run
start_time = time.perf_counter()



#useful mathematical function to produce array of prime factors of integer n but with any squares removed
def Remove_squares(n):
    factors = factorint(n)
    prime_factors = [] #i think normal array better for looping over than numpy array
    for factor, exponent in factors.items():
        if exponent % 2 == 1:  #only include factors with odd exponents
            prime_factors.append(factor) #outputs an array
        if prime_factors == []: #e.g. if the number is a square...
            prime_factors = [1]
    return np.prod(np.array(prime_factors)).astype(int)


#scan through the integers to find a solution. Not the most efficient or sophisticated method, 
#but it works fairly quickly for most simple quadratic forms.
def Solution_scanner(Q, max_value):
    n = Q.shape[0]
    for indices in np.ndindex(*([max_value] * n)):
        sol = np.array(indices) + 1  #add 1 to each index to get values between 1 and max_value
        if sol.T @ Q @ sol == 0:
            return sol
    raise ValueError("No solutions were found")



#Find two distinct solutions, thus defining a hyperbolic plane
def DistinctSolutions(Q, max_value, rand_max):
    n = Q.shape[0]
    #find v1 such that v1^T B v1 = 0
    v1 = Solution_scanner(Q, max_value)
    #find v2 such that v2^T B v2 = 0 and v2^T B v1 != 0 using the method of chords
    #start by generating a set of random integers we use as a parameterisation
    arr = np.random.randint(0, rand_max, size=n)
    phi = (v1-arr).T @ Q @ (v1-arr)
    psi = v1.T @ Q @ (arr-v1)
    #solution
    v2 = phi*v1 - 2*psi*(arr - v1)
    if v2.T @ Q @ v1 == 0:
        #this is unlikely. If it occurs, run again, and a new set of random variables will produce a new solution.
        raise ValueError("Try again. Not hyperbolic.")
        return None
    elif v2.T @ Q @ v2 == 0 and v1.T @ Q @ v1 == 0:
        return np.array([v1, v2])
        print("v1 is ", v1)
        print("v2 is ", v2)
    else:
        raise ValueError("Something has gone tragically wrong.")



#define function that generates orthogonal basis elements
def newfunc(v_basis, n):
    #n is number of elements to be generated
    #v_basis is the first 2 basis elements in the basis
    for i in range(2, n):
        #define the number of elements in the basis
        basis_size = len(v_basis)

        #define the matrix whose null space we want to find
        A = Matrix(np.vstack([v_basis[k].T @ Q for k in range(basis_size)]))
        
        #find the null space
        null_basis = A.nullspace()
        v = np.array(null_basis[0]).astype(np.float64).squeeze()
        
        #convert to integers
        denoms = [Fraction(x).limit_denominator().denominator if x is not None and not np.equal(np.mod(x, 1), 0) else 1 for x in v]

        lcm = np.lcm.reduce(denoms)
        v = (v * lcm).astype(int)
        
        #divide by the gcd
        nonzero_indices = np.nonzero(v)
        # Extract the non-zero elements
        nonzero_elements = v[nonzero_indices]
        gcd = np.gcd.reduce(nonzero_elements)
        #divide by gcd
        v = (v / gcd).astype(int)
        
        
        #add this vector to our basis vectors
        v_basis = np.vstack([v_basis, v])
        
    return v_basis
    





#define function to perform a Witt decomposition
def decompose(Q, max_value, rand_max):
    #find the dimension of Q
    n = Q.shape[0]
    
    #remove any common factors in Q
    diagonal = np.diag(Q)
    gcd = np.gcd.reduce(diagonal)
    
    Q = Q // gcd
    diagonal = np.diag(Q)
    dummy_vec = np.vectorize(Remove_squares)(diagonal)
    
    Q = np.diag(dummy_vec)
    
    print(Q)
    
    
    
    #find the Witt index. Don't forget to import witt_index.py!!
    s = Signature(Q)
    a = Absolute_determinant(Q)
    dimQ = Dim(Q)
    A = Adjusted_hasse_invariant(Q)
    d = Monster_function_d(s, a, A)
    W = Witt_index(dimQ, d)
    
    #run if we have isotropic vectors, and else return the anisotropic quadratic form
    if W == 0:
        return Q
    elif W > 0:
        print("Witt index is: ", W)
        #compute new basis
        #first two basis vectors are the given solutions
        v_basis = DistinctSolutions(Q, max_value, rand_max)
        v_basis = newfunc(v_basis, n)
        v_basis_T = v_basis.T
        
        #now re-express Q in this basis
        new_Q = v_basis @ Q @ v_basis_T

        return "xy ⊕ ", decompose(new_Q[2:, 2:], max_value, rand_max)
    else:
        raise ValueError("You fool!")
        
        
        
Q = np.diag([1,3,-2,-2])

print(decompose(Q, 10, 5))



"""runtime"""

#output time taken for code to run
end_time = time.perf_counter()
elapsed_time = end_time - start_time

print(f"Elapsed time: {elapsed_time:.4f} seconds")


