// #![feature(test)]

//! A lightweight Rust library for numerical integration of real-valued functions.
//!
//! `integrate` approximates definite integrals of the form
//!
//! $$\int_a^b f(x)\\, dx$$
//!
//! using Newton-Cotes formulas, Gaussian quadrature, Romberg's method, and
//! adaptive techniques.
//!
//! # Methods
//!
//! ## Newton-Cotes ([`newton_cotes`])
//!
//! | Function | Description |
//! |---|---|
//! | [`newton_cotes::rectangle_rule`] | Midpoint rectangle rule |
//! | [`newton_cotes::trapezoidal_rule`] | Trapezoidal rule |
//! | [`newton_cotes::simpson_rule`] | Simpson's 1/3 rule |
//! | [`newton_cotes::newton_rule`] | Newton's 3/8 rule |
//!
//! ## Gaussian Quadrature ([`gauss_quadrature`])
//!
//! | Function | Description |
//! |---|---|
//! | [`gauss_quadrature::legendre_rule`] | Gauss-Legendre rule |
//! | [`gauss_quadrature::gauss_laguerre_rule`] | Gauss-Laguerre rule |
//! | [`gauss_quadrature::gauss_hermite_rule`] | Gauss-Hermite rule |
//! | [`gauss_quadrature::gauss_first_kind_chebyshev_rule`] | Gauss-Chebyshev rule (first kind) |
//! | [`gauss_quadrature::gauss_second_kind_chebyshev_rule`] | Gauss-Chebyshev rule (second kind) |
//!
//! ## Adaptive Quadrature ([`adaptive_quadrature`])
//!
//! | Function | Description |
//! |---|---|
//! | [`adaptive_quadrature::adaptive_simpson_method`] | Adaptive Simpson's method |
//!
//! ## Romberg ([`romberg`])
//!
//! | Function | Description |
//! |---|---|
//! | [`romberg::romberg_method`] | Romberg integration |
//!
//! # Quick Start
//!
//! ```rust
//! use integrate::prelude::*;
//! let result = trapezoidal_rule(|x: f64| x.exp(), 0.0, 1.0, 1000_usize);
//! assert!((result - (std::f64::consts::E - 1.0)).abs() < 1e-6);
//! ```
//!
//! # Caveats
//!
//! All of the numerical integration techniques listed above assume that
//! the integrand is continuous on the interval of integration. The error
//! estimates generally require that the integrands are differentiable of
//! whatever order is required so that the formula for the error estimate
//! makes sense. Below is a checklist to verify before using any of these
//! algorithms.
//!
//! ### Integrand checklist
//!
//! - Are there any singularities of the integrand in the interval of integration?
//!
//!     - If there are, then can they be removed?
//!
//!     - A function which has jump discontinuities can be integrated by splitting
//!       the interval of integration into subintervals on which the function is
//!       continuous, then numerically integrating over each subinterval and summing
//!       the results.
//!
//!     - Using integration by parts, certain types of singular integrals can be
//!       expressed as the sum of a closed-form term and a non-singular integral.
//!
//!     - In other cases, an approximation of the integral in a small neighborhood of
//!       the singularity can be obtained analytically, and the interval of integration
//!       can be split so that numerical integration is used on the non-singular parts.
//!
//! - Does the function oscillate over the region of integration? If so, make sure
//!   that the step size is smaller than the wavelength of the function. The interval
//!   can also be split into subintervals of half a wavelength each.

pub mod adaptive_quadrature;
pub mod gauss_quadrature;
pub mod newton_cotes;
pub mod prelude;
pub mod romberg;
mod utils;
