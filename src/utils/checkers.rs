use num::{Float, Unsigned, Zero};

/// Validates arguments common to all Newton-Cotes integration methods.
///
/// * `a` - lower limit of the integration interval.
/// * `b` - upper limit of the integration interval.
/// * `n` - number of steps (sub-intervals).
///
/// # Panics
///
/// - Panics if `n` is zero.
/// - Panics if `a` or `b` is infinite or NaN.
/// - Panics if `a > b` (the lower limit must not exceed the upper limit).
pub fn check_newton_method_args<F: Float, U: Unsigned>(a: F, b: F, n: U) {
    if n.is_zero() {
        panic!("number of steps can't be zero");
    }

    if a.is_infinite() | b.is_infinite() {
        panic!("Integral limits a and b can't be infinite");
    }

    if a.is_nan() || b.is_nan() {
        panic!("Integral limits a and b can't be NaN");
    }

    if a > b {
        panic!("a must be strictly less than b");
    }
}

/// Validates the degree argument for Gaussian quadrature rules.
///
/// * `n` - number of quadrature nodes (must be at least 1).
///
/// # Panics
///
/// - Panics if `n` is zero.
pub fn check_gauss_rule_args(n: usize) {
    if n.is_zero() {
        panic!("number of steps can't be zero");
    }
}
