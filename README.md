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

The two methods are compared analytically in a Mathematica file for one example, and the sets of solutions of the two different parameterisations are compared
in a Python file for two examples. The overlap is relatively small in both cases, but this is unsurprising (there are, after all, infinitely many rational number with which to parameterise...).

I may upload the final report for this project at a later date.

## Technologies

* Python version: 3.9.14
* Mathematica version: 13.2

## List of Files
* **method-of-chords.py**: finds integer solutions to via parameterisation $x^2 + y^2 = 2z^2$
* **brute-force.py**: scans over integers for solutions to $x^2 + y^2 = 2z^2$
* **two-gauge-fields.py**:  scans for $n=2$ solutions to $x^2 + y^2 = 2z^2$
* **witt-index.py**: calculates the Witt index $W$ for an arbitrary rational quadratic form
* method-of-chords-2.py: an example of using the method of chords in a Witt decomposed system
        $$x^2 + y^2 - w^2 - 3z^2 = 0$$
* **orthogonal-matrices-4x4.py**: an example of using the generalised orthogonal matrix method in a Witt decomposed system $$x^2 + y^2 - w^2 - 3z^2 = 0$$ to generate generalised orthogonal matrices
* **orthogonal-matrices-final.nb**: parameterisation of orthogonal matrices to produce all maximal isotropic subspaces for examples with $W=1,2$
* **orthogonal-MoC-comparison.py**: comparison of solution sets produced by the Method of Chords and the generalised orthogonal matrix method
* **second-example-MoC-matrix-comparison.py**: comparison of solution sets in a second case
* **MoC-orthogonal-matrix-comparison-analytical.nb**: analytical comparison of the Method of Chords and the generalised orthogonal matrix method for $$x^2 + y^2 - w^2 - 3z^2 = 0$$ 



## References
<a id="1">[1]</a> 
Beale and Harrison (1989). 
A computation of the Witt index for rational quadratic forms.
Aequationes Mathematicae, Vol. 38, 86--98.
