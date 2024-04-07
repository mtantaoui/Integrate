//! Romberg's method is used to estimate the integral of a function on a closed
//! and bounded interval.
//!
//! Classically, the method consists of successively applying the composite
//! trapezoidal rule each time halving the length of the subintervals and
//! using a linear combination of the resulting sequence of estimates
//! to estimate the integral by successively deleting the low order error
//! terms in the Euler-Maclaurin summation formula.
//!
//! The process terminates when the change of the estimate is within a
//! preassigned tolerance within a preassigned number of successive estimates.

pub mod method;
