//! Romberg's method for numerical integration via Richardson extrapolation.
//!
//! Builds a triangular table of trapezoidal estimates at successively halved
//! step sizes, then applies Richardson extrapolation column by column to
//! accelerate convergence. The final entry $R[n-1, n-1]$ has error $O(h^{2n})$
//! for smooth integrands.

use num::{Float, ToPrimitive, Unsigned};

use crate::utils::checkers::check_newton_method_args;

/// Approximates $\int_a^b f(x)\\,dx$ using Romberg's method.
///
/// Romberg's method builds a triangular extrapolation table from trapezoidal
/// estimates at step sizes $h, h/2, h/4, \ldots$ The bottom-right entry
/// $R[n-1, n-1]$ gives the best estimate.
///
/// # Algorithm
///
/// **Phase 1 – first column.**  $R[0, 0]$ is the single-interval trapezoid.
/// Each subsequent row reuses the previous estimate and adds only the
/// $2^{i-1}$ new midpoints, so the total number of function evaluations is
/// $2^{n-1} + 1$ (not $O(n \cdot 2^n)$):
/// $$R[i, 0] = \tfrac{1}{2} R[i-1, 0] + h_i \sum_{k=0}^{2^{i-1}-1} f\!\left(a + (2k+1)h_i\right)$$
///
/// **Phase 2 – Richardson extrapolation.**  Fill left to right:
/// $$R[i, j] = \frac{4^j \, R[i, j-1] - R[i-1, j-1]}{4^j - 1}$$
///
/// Only two rows are kept in memory at any time.
///
/// # Parameters
/// * `func`        – integrand $f$, a function of a single variable.
/// * `lower_limit` – lower bound $a$ of the integration interval.
/// * `upper_limit` – upper bound $b$ of the integration interval.
/// * `n_columns`   – number of extrapolation columns; the method uses $2^{n-1}+1$
///   function evaluations and achieves $O(h^{2n})$ accuracy for smooth $f$.
///
/// # Panics
/// Panics if `n_columns` is zero, either limit is non-finite, or
/// `lower_limit > upper_limit`.
///
/// # Examples
/// ```
/// use integrate::romberg::romberg_method;
///
///
/// fn square(x: f64) -> f64 {
///     x.powi(2)
/// }
///
/// let a = 0.0;
/// let b = 1.0;
///
/// let num_steps: usize = 10;
///
/// let integral = romberg_method(square, a, b, num_steps);
/// ```
///
/// # Resources
/// * [Methods of numerical Integration (2nd edition), by Philip J. Davis and Philip Rabinowitz.](https://www.cambridge.org/core/journals/mathematical-gazette/article/abs/methods-of-numerical-integration-2nd-edition-by-philip-j-davis-and-philip-rabinowitz-pp-612-3650-1984-isbn-0122063600-academic-press/C331158D0392E1D5CD9B0C6ED4EE5F43)
/// * [Romberg's method](https://en.wikipedia.org/wiki/Romberg%27s_method)
pub fn romberg_method<Func, F1, F2, U>(
    func: Func,
    lower_limit: F1,
    upper_limit: F1,
    n_columns: U,
) -> f64
where
    F1: Float,
    F2: Float,
    U: Unsigned + ToPrimitive + Copy,
    Func: Fn(F1) -> F2,
{
    check_newton_method_args(lower_limit, upper_limit, n_columns);

    let n = n_columns.to_usize().unwrap();
    let a = lower_limit.to_f64().unwrap();
    let b = upper_limit.to_f64().unwrap();
    let h0 = b - a;

    // Wrap func to work in f64.
    let f = |x: f64| func(F1::from(x).unwrap()).to_f64().unwrap();

    // Rolling two-row buffer: prev = row i-1, curr = row i.
    let mut prev = vec![0.0_f64; n];
    let mut curr = vec![0.0_f64; n];

    // R[0, 0]: single-interval trapezoidal estimate.
    prev[0] = 0.5 * h0 * (f(a) + f(b));

    for i in 1..n {
        let num_new = 1_usize << (i - 1); // 2^(i-1) new midpoints
        let h = h0 / (2 * num_new) as f64; // step size at level i

        // R[i, 0]: reuse prev[0] and add only the new midpoints.
        let new_sum: f64 = (0..num_new).map(|k| f(a + (2 * k + 1) as f64 * h)).sum();
        curr[0] = 0.5 * prev[0] + h * new_sum;

        // Richardson extrapolation: fill columns 1..=i.
        let mut factor = 4.0_f64; // 4^j, starting at j=1
        for j in 1..=i {
            curr[j] = (factor * curr[j - 1] - prev[j - 1]) / (factor - 1.0);
            factor *= 4.0;
        }

        std::mem::swap(&mut prev, &mut curr);
    }

    prev[n - 1]
}

#[cfg(test)]
mod tests {
    use std::ops::Div;

    use super::*;

    const EPSILON: f64 = 10e-5;
    const NUM_STEPS: usize = 10;

    #[test]
    fn test_integral_value() {
        fn square(x: f64) -> f64 {
            x.powi(2)
        }

        let a = 0.0;
        let b = 1.0;

        let integral = romberg_method(square, a, b, NUM_STEPS);

        let analytic_result: f64 = 1.0.div(3.0);

        assert!((integral - analytic_result).abs() < EPSILON);
    }

    #[test]
    fn test_f32_to_f64() {
        fn square(x: f32) -> f64 {
            x.powi(2) as f64
        }

        let a = 0.0;
        let b = 1.0;

        let integral = romberg_method(square, a, b, NUM_STEPS);

        let analytic_result: f64 = 1.0.div(3.0);

        assert!((integral - analytic_result).abs() < EPSILON);
    }

    #[test]
    fn test_f64_to_f32() {
        fn square(x: f64) -> f32 {
            x.powi(2) as f32
        }

        let a = 0.0;
        let b = 1.0;

        let integral = romberg_method(square, a, b, NUM_STEPS);

        let analytic_result: f64 = 1.0.div(3.0);

        assert!((integral - analytic_result).abs() < EPSILON);
    }

    #[test]
    fn test_f32_to_f32() {
        fn square(x: f32) -> f32 {
            x.powi(2)
        }

        let a = 0.0;
        let b = 1.0;

        let integral = romberg_method(square, a, b, NUM_STEPS);

        let analytic_result: f64 = 1.0.div(3.0);

        assert!((integral - analytic_result).abs() < EPSILON);
    }

    #[test]
    #[should_panic(expected = "Integral limits a and b can't be NaN")]
    fn test_lower_limit_nan() {
        fn f(x: f64) -> f64 {
            x
        }
        let _ = romberg_method(f, f64::NAN, 1.0, 1usize);
    }

    #[test]
    #[should_panic(expected = "Integral limits a and b can't be NaN")]
    fn test_upper_limit_nan() {
        fn f(x: f64) -> f64 {
            x
        }
        let _ = romberg_method(f, 0.0, f64::NAN, 1usize);
    }
}
