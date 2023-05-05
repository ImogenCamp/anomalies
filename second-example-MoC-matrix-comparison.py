#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  3 17:10:02 2023

@author: imogencamp
"""
"""This code does another comparison between the method of chords and the generalised orthogonal matrix method.
This code was written to check the results in the previous project (Michael's). It probably also uses a better
method of removing duplicate solutions...

Note the overlap here is much smaller (approx. 1%).
"""


import numpy as np
import time 

#we wish to time how long it takes the code to run
start_time = time.perf_counter()

#function to check if diiophantine equation is satisfied
def Diophantine_checker(B, sol):
    squared_sol = np.power(sol, 2)
    diophantine_sum = np.dot(B, squared_sol)
    if diophantine_sum == 0:
        return True
    else:
        return [False]
        

"""orthogonal matrix method"""

min_value = 1
max_value = 10

B = np.array([1, 1, -1, -1])
sol_array = []
truth_array = []

for a in range(min_value, max_value):
    for b in range(min_value, max_value):
        for c in range(min_value, max_value):
            for d in range(min_value, max_value):
                for e in range(min_value, max_value):
                    for f in range(min_value, max_value):
                            solution = np.array([1 - 2*c + c**2 - d**2 + 2*c*d**2 - c**2*d**2 - e**2 - b**2*(-1 + e**2) + 
                                                 f**2 + 2*b*((-1 + c)*d*e + f) + 2*a*(e + (-1 + c)*d*f - b*e*f) - a**2*(1 + f**2), 
                                                 
                                                 2*(b**2*e + (-1 + c)*(e - d*f) + b*(d - c*d + e*f) 
                                                 + a*(1 - c + b*f + f**2)), 
                                                 
                                                 -2*(-b*(-1 + c - a*e + e**2) + (-1 + c)*(d*e - f) 
                                                + a**2*f + a*(d - c*d - e*f)), 
                                                                                                           
                                                1 - 2*c + c**2 - d**2 + 2*c*d**2 - c**2*d**2 + e**2 - b**2*(1 + e**2) 
                                                - f**2 - 2*b*(d*e - c*d*e + f) - 2*a*(e - (-1 + c)*d*f + b*e*f)
                                                - a**2*(-1 + f**2)])
                            
                            gcd_val = np.gcd.reduce(np.abs(solution))
                            
                            if gcd_val != 0:
                                solution = solution // gcd_val
                            
                            sol_array.append(solution)
                            truth_array.append(Diophantine_checker(B, solution))
                                                                                 

                            
sol_array = np.array(sol_array)
truth_array = np.array(truth_array)

#print solutions and check they are solutions
'''
if np.all(truth_array):
    print("All elements are True")
    print(sol_array.shape[0])
    print(sol_array)
else:
    print("Not all elements are True")
'''
    
"""method of chords"""

min_val = 1
max_val = 100

moc_sol = []
moc_truth = []

for a in range(min_val, max_val):
    for b in range(min_val, max_val):
        for c in range(min_val, max_val):
            solution = np.array([-a**2 + b**2 - c**2, -2*a*b, -2*a*c, a**2 + b**2 - c**2])
            
            gcd_val = np.gcd.reduce(np.abs(solution))
            solution = solution // gcd_val
            
            
            moc_sol.append(solution)
            moc_truth.append(Diophantine_checker(B, solution))
            

moc_sol = np.array(moc_sol)
moc_truth = np.array(moc_truth)

#print solutions and check they are solutions
'''
if np.all(moc_truth):
    print("All elements are True")
    print(moc_sol.shape[0])
    print(moc_sol)
else:
    print("Not all elements are True")
'''

moc_tuples = tuple(map(tuple, moc_sol))
orth_tuples = tuple(map(tuple, sol_array))

moc_set = set(tup for tup in moc_tuples)
orth_set = set(tup for tup in orth_tuples)

print("number of MoC solutions: ", len(moc_set))
print("number of orthogonal matrix method solutions: ", len(orth_set))

#print("solutions from orthogonal matrix method: ", orth_set)
#print("solutions from MoC: ", moc_set)

common_solutions = orth_set.intersection(moc_set)
print("number of common solutions: ", len(common_solutions))
#print(common_solutions)
#print(orth_set - common_solutions)



"""runtime"""

#output time taken for code to run
end_time = time.perf_counter()
elapsed_time = end_time - start_time

print(f"Elapsed time: {elapsed_time:.4f} seconds")


"""Output:
    
    number of MoC solutions:  811213
    number of orthogonal matrix method solutions:  117015
    number of common solutions:  1547
    Elapsed time: 17.3901 seconds

"""













   