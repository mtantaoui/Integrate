# Integrator

a small, lightweight numerical integration written in Rust.

# TODO:

- gershgorin_bound method switched does not check that first element of subdiagonal is zero, how to solve this issue ?
- implement https://www.netlib.org/lapack/lawnspdf/lawn132.pdf
- better error handling for matrices.
- write tests for matrices.
- write documentation for matrices.
- compare sequential strum sequence to parallelized strun sequence.
- parallelize givens bisection algorithm

# references for Laguerre polynomial zeros

- http://naturalunits.blogspot.com/2014/06/zeros-of-laguerre-polynomial.html
- https://en.wikipedia.org/wiki/Laguerre_polynomials#Recursive_definition,_closed_form,_and_generating_function (Recursive definition, closed form, and generating function)
- compare (approximate + ITP method) of zero-finding to computing eigenvalues of companion matrix method.

# references for Givens bisection algorithm:

- https://math.mit.edu/~edelman/homepage/papers/sturm.pdf
- https://ia902301.us.archive.org/2/items/c-36_20211010/C36.pdf
