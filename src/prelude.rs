//! The integrate prelude imports the various Numerical Integration methods.
//! The intention is that one can include `use integrate::prelude::*` and
//! have easy access to the various methods you will need.

pub use crate::adaptive_quadrature::adaptive_simpson_method;
pub use crate::gauss_quadrature::{
    gauss_first_kind_chebyshev_rule, gauss_hermite_rule, gauss_laguerre_rule,
    gauss_second_kind_chebyshev_rule, legendre_rule,
};
pub use crate::newton_cotes::{newton_rule, rectangle_rule, simpson_rule, trapezoidal_rule};
pub use crate::romberg::romberg_method;
