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

pub mod rectangle;
