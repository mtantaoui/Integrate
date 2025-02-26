use num::{Float, Unsigned, Zero};

/// Checks integral arguments for Newton-Codes methods
///
/// * `a` - lower limit of the integration interval.
/// * `b` - lower limit of the integration interval.
/// * `n` - number of steps.
pub fn check_newton_method_args<F: Float, U: Unsigned>(a: F, b: F, n: U) {
    if n.is_zero() {
        panic!("number of steps can't be zero");
    }

    if a.is_infinite() | b.is_infinite() {
        panic!("Integral limits a and b can't be infinite");
    }

    if a > b {
        panic!("a must be strictly less than b");
    }
}

/// Checks integral arguments for Gauss-Laguerre rule
///
/// * `n` - number of steps.
pub fn check_gauss_rule_args(n: usize) {
    if n.is_zero() {
        panic!("number of steps can't be zero");
    }
}
