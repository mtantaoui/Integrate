//! Newton-Cotes methods
//!
//! Newton-Cotes methods approximate the integral of a function by summing a weighted
//! combination of the function evaluated at equally-spaced points and nodes. If the
//! endpoints of the interval of integration are excluded from the sum, the method
//! is called an open Newton-Cotes method and if the endpoints are included the
//! method is called a closed Newton-Cotes method.
//!
//! - Rectangle Rule.
//! - Trapezoidal Rule.
//! - Simpson's Rule.
//! - Newton's 3/8 Rule.

pub mod newton;
pub mod rectangle;
pub mod simpson;
pub mod trapezoidal;
mod utils;
