#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 19:22:40 2023

@author: imogencamp
"""

"""The goal of this code is to implement the algorithm given in Harrison & Beale, allowing us to calculate
the Witt index of a rational quadratic form. This is used on some examples.
"""

import numpy as np
import math
from sympy.ntheory import legendre_symbol
from sympy import isprime
from sympy import factorint

#we have a diagonal quadratic form f(x)
#can be expressed as x^T Q x for a diagonal matrix Q with no zero elements
#we will consider Q to be our input
#algorithm follows from paper beale & harrison


"""generally useful mathematical functions"""

#function which takes symmetric difference of two sets exists already: A.symmetric_difference(B)
#print({3,4}.symmetric_difference({3}))

#useful mathematical function to produce array of prime factors of integer n but with any squares removed
def Prime_factors_no_squares(n):
    factors = factorint(n)
    prime_factors = [] #i think normal array better for looping over than numpy array
    for factor, exponent in factors.items():
        if exponent % 2 == 1:  #only include factors with odd exponents
            prime_factors.append(factor) #outputs an array
        if prime_factors == []: #e.g. if the number is a square...
            prime_factors = [1]
    return np.array(prime_factors)


#legendre symbol exists in Python but we must modify to make consistent with beale & harrison
def Legendre_symbol(a, prime):
    l = legendre_symbol(a, prime)
    if l == -1:
        return 1
    else:
        return 0
    
#print(legendre_symbol(2, 3))
#print(Legendre_symbol(2, 3))


"""functions for algorithm"""

#witt index
def Witt_index(dimQ, d):
    return round((dimQ - d)/2)


#signature of the quadratic form
def Signature(Q):
    return np.count_nonzero(Q > 0) - np.count_nonzero(Q < 0)


#dimension of the quadratic form
def Dim(Q):
    return np.count_nonzero(Q > 0) + np.count_nonzero(Q < 0)


#absolute determinant of quadratic form
def Absolute_determinant(Q):
    # Find the determinant of the matrix
    det = round(abs(np.linalg.det(Q)))
    
    #remove the square factors in the determinant
    for i in range(2, int(np.sqrt(det)) + 1):
        if det % (i**2) == 0:
            det = det // (i**2)
            
    return det



#define symmetric bilinear map
#technically this is defined over the rationals, but our inputs will only ever be positive integers...
#hence, we will be lazy <3
def Phi(p, q):
    if isprime(p) and isprime(q):
        if p == 2 and q == 2:
            return set() #empty set
        elif p != 2 and q == 2:
            n = Legendre_symbol(2, p)
            return {p} if n != 0 else set()
        elif p == 2 and q != 2:
            n = Legendre_symbol(2, q)
            return {q} if n != 0 else set()
        elif p == q:
            n = Legendre_symbol(-1, p)
            return {p} if n != 0 else set()
        elif p != q:
            n = Legendre_symbol(q, p)
            m = Legendre_symbol(p, q)
            set1 = {p} if n != 0 else set()
            set2 = {q} if m != 0 else set()
            return set1.symmetric_difference(set2)
    elif (p == 1) or (q == 1):
        return set() #empty set
    elif (not isprime(p)) or (not isprime(q)):
        #we remove the squares since these will cancel when taking the symmetric difference
        factors_p = Prime_factors_no_squares(p)
        factors_q = Prime_factors_no_squares(q)
        initial_set = set()
        #use properties of symmetric biinear map to show step below makes sense
        for a in factors_p:
            for b in factors_q:
                initial_set = initial_set.symmetric_difference(Phi(a, b))
        return initial_set
    else:
        #something will have gone wrong, print error message
        print("ERROR IN PHI")
        
#print(Phi(10,2))
        

#adjusted Hasse invariant of the quadratic form
def Adjusted_hasse_invariant(Q):
    diag_elements = np.diag(Q)
    dimension = Dim(Q)
    first_term = set()
    second_term = set()
    
    # Calculate the first term
    pairs = [(abs(diag_elements[i]), abs(diag_elements[j])) for i in range(dimension) for j in range(i+1, dimension)]
    for pair in pairs:
        first_term = first_term.symmetric_difference(Phi(pair[0], pair[1]))
    
    # Calculate the second term
    for j in range(dimension):
        if diag_elements[j] > 0:
            second_term = second_term.symmetric_difference(Phi(diag_elements[j], diag_elements[j]))

    return second_term.symmetric_difference(first_term)
        
    
#define the really nasty function from beale & harrison        
#not too much commenting, since not too much to add...
def Monster_function_d(s, a, A):
    if abs(s) > 2:
        return abs(s)
    elif s == 2:
        if any(math.gcd(p, a) == 1 and Legendre_symbol(-a, p) == 1 for p in A):
            return 4
        elif any(p % 8 == 7 for p in A) and len(A) % 2 == 0:
            return 4
        else: 
            return 2
    elif s == 1:
        if Phi(a, a) == A:
            return 1
        else:
            return 3
    elif s == 0 and a == 1 and A == set():
        return 0
    elif s == 0 and not (a == 1 and A == set()):
        if any(math.gcd(p, a) == 1 and Legendre_symbol(a, p) == 1 for p in A):
            return 4
        elif any(p % 8 == 1 for p in A) and len(A) % 2 == 1:
            return 4
        else: 
            return 2
    elif s == -1:
        if A == set():
            return 1
        else:
            return 3
    elif s == -2:
        if any(math.gcd(p, a) == 1 and Legendre_symbol(-a, p) == 1 for p in A):
            return 4
        elif any(p % 8 == 7 for p in A) and len(A) % 2 == 0:
            return 4
        else: 
            return 2
    else:
        #something will have gone wrong, print an error message
        print("ERROR IN MONSTER FUNCTION")

"""Now check this is all working with some examples"""

Q1 = np.diag([1, 1, -2])

s1 = Signature(Q1)
a1 = Absolute_determinant(Q1)
dimQ1 = Dim(Q1)
A1 = Adjusted_hasse_invariant(Q1)
d1 = Monster_function_d(s1, a1, A1)

#witt index should be 1
w1 = Witt_index(dimQ1, d1)
print(s1, a1, dimQ1, A1, d1, "Witt index is:", w1)
#YAY IT WORKS

Q2 = np.diag([1, 1, -1, -1])

s2 = Signature(Q2)
a2 = Absolute_determinant(Q2)
dimQ2 = Dim(Q2)
A2 = Adjusted_hasse_invariant(Q2)
d2 = Monster_function_d(s2, a2, A2)

#witt index should be 2
w2 = Witt_index(dimQ2, d2)
print(s2, a2, dimQ2, A2, d2, "Witt index is:", w2)
#YAY IT WORKS


Q3 = np.diag([3, 7, -4, -5, -1])

s3 = Signature(Q3)
a3 = Absolute_determinant(Q3)
dimQ3 = Dim(Q3)
A3 = Adjusted_hasse_invariant(Q3)
d3 = Monster_function_d(s3, a3, A3)

#initial scan would suggest witt index of 1
w3 = Witt_index(dimQ3, d3)
print(s3, a3, dimQ3, A3, d3, "Witt index is:", w3)
#this is indeed the case



Q3 = np.diag([1, 4, 5,-1,-7,-2])

s3 = Signature(Q3)
a3 = Absolute_determinant(Q3)
dimQ3 = Dim(Q3)
A3 = Adjusted_hasse_invariant(Q3)
d3 = Monster_function_d(s3, a3, A3)


w3 = Witt_index(dimQ3, d3)
print(s3, a3, dimQ3, A3, d3, "Witt index is:", w3)



Q3 = np.diag([1, 1, 6])

s3 = Signature(Q3)
a3 = Absolute_determinant(Q3)
dimQ3 = Dim(Q3)
A3 = Adjusted_hasse_invariant(Q3)
d3 = Monster_function_d(s3, a3, A3)


w3 = Witt_index(dimQ3, d3)
print(s3, a3, dimQ3, A3, d3, "Witt index is:", w3)



Q3 = np.diag([1,3,-2,-2])

s3 = Signature(Q3)
a3 = Absolute_determinant(Q3)
dimQ3 = Dim(Q3)
A3 = Adjusted_hasse_invariant(Q3)
d3 = Monster_function_d(s3, a3, A3)


w3 = Witt_index(dimQ3, d3)
print(s3, a3, dimQ3, A3, d3, "Witt index is:", w3)








