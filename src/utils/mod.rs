//! Internal utilities for polynomial root-finding, quadrature node/weight
//! computation, and input validation.
//!
//! This module is **not** part of the public API.

pub mod adaptive_simpson;
pub mod bessel;
pub mod chebyshev;
pub mod checkers;
pub mod hermite;
pub mod laguerre;
pub mod legendre;
pub mod matrix;
pub mod orthogonal_polynomials;
