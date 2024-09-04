use num::Zero;

/// Checks integral arguments for Gauss-Laguerre rule
///
/// * `n` - number of steps.
pub fn check_gauss_rule_args(n: usize) {
    if n.is_zero() {
        panic!("number of steps can't be zero");
    }
}
