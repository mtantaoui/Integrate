# Adaptive Simpson

If an integrand is poorly behaved in a small interval about a point,
then an attempt to integrate the function over an interval which contains
the poorly behaved interval either requires that small subintervals
are chosen for composite quadratures or the interval is decomposed into three intervals,
two on which the function is well-behaved and relatively large subintervals
can be chosen for the composite quadrature technique and one in which smaller subintervals need to be chosen.

Adaptive techniques are attempts to automatically detect and control the length of subintervals.

The technique for which the link to the listing is given below uses Simpson's rule
for integrating a function \\(f(x)\\) on a closed and bounded interval \\(\[a,b\]\\).
