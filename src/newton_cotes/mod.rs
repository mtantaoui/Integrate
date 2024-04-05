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

use num::{Float, Unsigned};

pub mod newton;
pub mod rectangle;
pub mod simpson;
pub mod trapezoidal;

/// Checks integral arguments for Newton-Codes methods
pub fn check_integral_args<F: Float, U: Unsigned>(a: F, b: F, n: U) {
    if n.is_zero() {
        panic!("number of steps can't be zero");
    }

    if a.is_infinite() | b.is_infinite() {
        panic!("Integral limits a and b can't be infinite");
    }
}
