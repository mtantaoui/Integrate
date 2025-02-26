# Gaussian quadrature

Instead of evaluating the integrand at equally spaced nodes as in Newton-Cotes methods,
Gaussian quadrature methods make a judicious choice of nodes so as to maximize the precision
of the numerical integration relative to the number of integrand evaluations.

The common Gaussian quadrature methods are:

- [Gauss-Legendre](gauss_quadrature/gauss_legendre.md) used to integrate a function \\(f(x)\\) over a closed and bounded interval \\(\[a,b\]\\).
- [Gauss-Laguerre](gauss_quadrature/gauss_laguerre.md) used to integrate a function of the form \\(f(x) e^{-x}\\) over the positive x-axis \\(\lbrace x \in \mathbb{R} : x > 0 \rbrace\\).
- [Gauss-Hermite](gauss_quadrature/gauss_hermite.md) used to integrate a function of the form \\(f(x) e^{-x^2}\\) over the entire x-axis, \\(\lbrace x \in \mathbb{R} : -\infty < x < \infty \rbrace\\).
- [Gauss-Chebyshev](gauss_quadrature/gauss_chebyshev.md) First Kind used to integrate a function of the form \\(\frac{f(x)}{\sqrt( 1-x^2 )}\\) over the interval \\(\[-1,1\]\\).
- [Gauss-Chebyshev](gauss_quadrature/gauss_chebyshev.md) Second Kind used to integrate a function of the form \\(f(x) \* \sqrt{ 1-x^2 }\\) over the interval \\(\[-1,1\]\\).
