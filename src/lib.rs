//! Numerical integration is concerned with developing algorithms to
//! approximate the integral of a function $f(x)$. The most commonly used algorithms
//! are Newton-Cotes formulas, Romberg's method, Gaussian quadrature, and to
//! lesser extents Hermite's formulas and certain adaptive techniques.
//!
//! # Caveats
//!
//! All of the numerical integration techniques listed above assume that
//! the integrand is continuous on the interval of integration. The error
//! estimates generally require that the integrands are differentiable of
//! whatever order is required so that the formula for the error estimate
//! makes sense. Below is a check list which one should verify before using
//! any of the numerical algorithms.
//!
//! ### Integrand check list
//!
//! - Are there any singularities of the integrand in the interval of integration?
//!
//!     - If there are, then can they be removed?
//!
//!     - A function which has jump discontinuities can be integrated by splitting
//!  the interval of integration into subintervals on which the function is continuous,
//!  and then numerically integrate over each subinterval and add the results.
//!     
//!     - Using integration by parts certain types of singular integrals can be
//! expressed as the sum of a function and a non-singular integral.
//!     
//!     - In other cases, an approximation of the integral in a small neighborhood of
//! the singularity can be obtained by hand and the interval of integration can be
//! split into three intervals in which numerical integration can be used on
//! two of them.
//!
//! - Does the function oscillate over the region of integration? If it does,
//!  then make sure that the step size is chosen to be smaller than the wave length
//! of the function. The interval of integration can also be split into subintervals
//! in which each subinterval is half a wave length and the algorithm is applied
//! to each subinterval.
//!
//! # Newton-Cotes Methods
//!
//! Newton-Cotes methods approximate the integral of a function by summing a weighted
//! combination of the function evaluated at equally-spaced points, nodes. If the
//! endpoints of the interval of integration are excluded from the sum, the method
//! is called an open Newton-Cotes method and if the endpoints are included the
//! method is called a closed Newton-Cotes method.
//!
//! - [x] Rectangle Rule.
//! - [] Trapezoidal Rule.
//! - [] Simpson's Rule.
//! - [] Newton's 3/8 Rule.

pub mod newton_cotes;