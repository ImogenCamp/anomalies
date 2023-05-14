# Anomaly cancellation in 2D quantum field theories

## Table of contents
* [General info](#general-info)
* [Technologies](#technologies)
* [List of Files](#list-of-files)
* [References](#references)


## General info

This is supplementary code written for this MSci project. The anomaly cancellation equations for a general non-Abelian group $G = H \times U(1)^n $
in 2D are systems of homogenous quadratic Diophantine equations for the $n$ sets of charges.
A single Diophantine equation ($n=1$) can be solved using the Method of Chords. Code was produced to generate solutions using this method. Code was also produced
using a brute force method to compare efficiency. Unsurprisingly, the brute force scan was hugely less efficient.

For larger systems of equations, solutions may or may not exist. Code was produced to scan through solutions and check if any linearly independent pair
satisfies the $n=2$ case. None were found.

The Theory of Quadratic Forms gives us that the maximal totally isotropic subspace has a dimension equal to the Witt index $W$. This is the maximum number
of linearly independent charges that can exist. An algorithm to calculate the Witt index of a rational quadratic form was produced in [[1]](#1). 
This was implemented to find the maximum number of linearly independent charges.

The Theory of Quadratic Forms brings to light another method for calculating solutions: acting on our initial solutions with elements of a generalisation
of the orthogonal group $$O = \\{\mathsf{O} \in \mathcal{M}(\mathbb{Q}) | \mathsf{O}^T \mathsf{B} \mathsf{O} = \mathsf{B}\\}.$$ These elements can be generated
using a generalisation of the Cayley parameterisation, which is done in a Mathematica file for two examples. Note that this holds for any $n$, provided we
have an initial solution of the correct dimension.

To generate such an initial solution, an algorithm was developed to explicitly perform a Witt decomposition. This is loosely based upon [[2]](#2), but is modified such that we have 2 initial solutions, rather than 1. This allows us to extract hyperbolic planes more easily. The code also does not use such sophisticated methods for generating solutions, as a simple scan suffices for our cases of interest.

The two methods are compared analytically in a Mathematica file for one example, and the sets of solutions of the two different parameterisations are compared
in a Python file for two examples. The overlap is relatively small in both cases, but this is unsurprising (there are, after all, infinitely many rational numbers with which to parameterise...).


## Technologies

* Python version: 3.9.14
* Mathematica version: 13.2

## List of Files
* **method-of-chords.py**: Finds integer solutions to via parameterisation $x^2 + y^2 = 2z^2$.
* **brute-force.py**: Scans over integers for solutions to $x^2 + y^2 = 2z^2$.
* **two-gauge-fields.py**:  Scans for $n=2$ solutions to $x^2 + y^2 = 2z^2$.
* **witt-index.py**: Calculates the Witt index $W$ for an arbitrary rational quadratic form.
* **method-of-chords-2.py**: An example of using the method of chords in a Witt decomposed system
        $$x^2 + y^2 - w^2 - 3z^2 = 0.$$
* **maximal_isotropic_subspace.py**: Implementation of an algorithm developed to explicitly perform a Witt decomposition.
* **orthogonal-matrices-4x4.py**: An example of using the generalised orthogonal matrix method in a Witt decomposed system $$x^2 + y^2 - w^2 - 3z^2 = 0$$ to generate generalised orthogonal matrices.
* **orthogonal-matrices-final.nb**: Parameterisation of orthogonal matrices to produce all maximal isotropic subspaces for examples with $W=1,2$.
* **orthogonal-MoC-comparison.py**: Comparison of solution sets produced by the Method of Chords and the generalised orthogonal matrix method.
* **second-example-MoC-matrix-comparison.py**: Comparison of solution sets in a second case.
* **MoC-orthogonal-matrix-comparison-analytical.nb**: Analytical comparison of the Method of Chords and the generalised orthogonal matrix method for $$x^2 + y^2 - w^2 - 3z^2 = 0.$$



## References
<a id="1">[1]</a> 
Beale and Harrison (1989). 
A computation of the Witt index for rational quadratic forms.
*Aequationes Mathematicae*, Vol. 38, 86--98.

<a id="2">[2]</a> 
D. Simon (2005). 
Quadratic equations in dimensions 4, 5 and more.
*Preprint*.
