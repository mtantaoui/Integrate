//! Convenience re-exports of all public integration functions.
//!
//! Import everything at once with:
//!
//! ```rust
//! use integrate::prelude::*;
//! ```
//!
//! # Re-exported symbols
//!
//! | Symbol | Description |
//! |---|---|
//! | [`adaptive_simpson_method`] | Adaptive Simpson's method with automatic step-size refinement |
//! | [`gauss_first_kind_chebyshev_rule`] | Gauss-Chebyshev quadrature rule of the first kind |
//! | [`gauss_hermite_rule`] | Gauss-Hermite quadrature rule for integrals over $(-\infty, +\infty)$ |
//! | [`gauss_laguerre_rule`] | Gauss-Laguerre quadrature rule for integrals over $[0, +\infty)$ |
//! | [`gauss_second_kind_chebyshev_rule`] | Gauss-Chebyshev quadrature rule of the second kind |
//! | [`legendre_rule`] | Gauss-Legendre quadrature rule |
//! | [`newton_rule`] | Newton's 3/8 rule |
//! | [`rectangle_rule`] | Midpoint rectangle rule |
//! | [`simpson_rule`] | Simpson's 1/3 rule |
//! | [`trapezoidal_rule`] | Trapezoidal rule |
//! | [`romberg_method`] | Romberg integration using Richardson extrapolation |

pub use crate::adaptive_quadrature::adaptive_simpson_method;
pub use crate::gauss_quadrature::{
    gauss_first_kind_chebyshev_rule, gauss_hermite_rule, gauss_laguerre_rule,
    gauss_second_kind_chebyshev_rule, legendre_rule,
};
pub use crate::newton_cotes::{newton_rule, rectangle_rule, simpson_rule, trapezoidal_rule};
pub use crate::romberg::romberg_method;
