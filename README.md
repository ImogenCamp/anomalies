# Anomaly cancellation in 2D quantum field theories

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
in a Python file for two examples. The overlap is relatively small in both cases, but this is unsurprising.

I may upload the final report for this project at a later date.






## References
<a id="1">[1]</a> 
Beale and Harrison (1989). 
A computation of the Witt index for rational quadratic forms.
Aequationes Mathematicae, Vol. 38, 86--98.
