//! Instead of evaluating the integrand at equally spaced nodes as in Newton-Cotes methods,
//! Gaussian quadrature methods make a judicious choice of nodes so as to maximize the precision
//! of the numerical integration relative to the number of integrand evaluations.
//! The common Gaussian quadrature methods are: Gauss-Legendre used to integrate a function $f(x)$
//! over a closed and bounded interval $\[a,b\]$, Gauss-Laguerre used to integrate a function of
//! the form $f(x) e^{-x}$ over the positive x-axis $\lbrace x \in \mathbb{R} : x > 0 \rbrace$,
//! Gauss-Hermite used to integrate a function of the form $f(x) e^{-x^2}$ over the entire x-axis, $\lbrace x \in \mathbb{R} : -\infty < x < \infty \rbrace$,
//! and Gauss-Chebyshev used to integrate a function of the form $\frac{f(x)}{\sqrt( 1-x^2 )}$ over the interval $\[-1,1\]$.
//!

mod bessel;
pub mod laguerre;
pub mod legendre;
