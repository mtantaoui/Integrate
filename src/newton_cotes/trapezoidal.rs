//! Trapezoidal Rule
//!
//! The trapezoidal rule approximates the integral of a function $f(x)$ on the closed and
//! bounded interval $\[a, a+h\]$ of length $h > 0$ by the (signed) area of the trapezoid formed
//! by the line segments joining $(a, 0)$ to $(a+h, 0)$, $(a+h, 0)$ to $(a+h, f(a+h))$, $(a+h, f(a+h))$
//! to $(a, f(a))$ and $(a, f(a))$ to $(a, 0)$.
//!
//! The composite trapezoidal rule is used to approximate the integral of a function
//! $f(x)$ over a closed and bounded interval $\[a, b\]$ where $a < b$,
//! by decomposing the interval $\[a, b\]$ into $n > 1$ subintervals of equal length
//! $h = \dfrac{b-a}{n}$, then adding the results of applying the trapezoidal
//! rule to each subinterval.
//!
//! By abuse of language both the composite trapezoidal rule and the trapezoidal rule
//! sometimes are referred to simply as the trapezoidal rule.
//!
//! Let $\int_{a}^{b} f(x) dx$ be the integral of $f(x)$ over the closed and bounded
//! interval $\[a,b\]$, and let $T_h(f)$ be the result of applying the trapezoidal
//! rule with $n$ subintervals of length h, i.e.
//!
//! ```math
//! T_h(f) = h \left[ \frac{f(a)}{2} + f(a+h) + ··· + f(b-h) + \frac{f(b)}{2} \right]
//! ```
//!
//! The Euler-Maclaurin summation formula relates $\int_{a}^{b} f(x) dx$ and $T_h(f)$
//! ```math
//! \begin{split}
//! T_h(f) &= \int_{a}^{b} f(x) dx + \frac{h^2}{12} \left[f'(b) - f'(a)\right] - \frac{h^4}{720} \left[f^{(3)}(b) - f^{(3)}(a) \right] \\
//! &+ ... + K h^{2p-2} \left[f^{(2p-3)}(b) - f^{(2p-3)(a)}\right] + O(h^{2p})
//! \end{split}
//! ```
//! where $f^\prime$, $f^{(3)}$, and $f^{(2p-3)}$ are the first, third and $(p-3)^{th}$ derivatives
//! of $f$ and $K$ is a constant.
//!
//! The last term, $O(h^{2p})$ is important. Given an infinitely differentiable function
//! in which the first $(2p-3)$ derivatives vanish at both endpoints of the interval of integration,
//! it is not true that $T_h(f) = \int_{a}^{b} f(x) dx$, but rather what the theorem says is that
//! ```math
//! \lim_{h \to 0} \left| \frac{T_h(f) - \int_{a}^{b} f(x)dx}{h^{2p}} \right| < M
//! ```
//! where $M>0$.
//!
//! If $f$ is at least twice differentiable on the interval $\[a,b\]$, then applying the mean-value theorem to
//! ```math
//! \begin{split}
//! T_h(f) - \int_{a}^{b} f(x) dx &= \frac{h^2}{12} \left[f^\prime(b) - f^\prime(a) \right] - \frac{h^4}{720} \left[ f^{(3)}(b) - f^{(3)}(a) \right] \\
//! &+ ... + K h^{2p-2} \left[f^{(2p-3)}(b) - f^{(2p-3)(a)} \right] + O(h^{2p})
//! \end{split}
//! ```
//! yields the standard truncation error expression
//!
//! ```math
//! R_h(f) - \int_{a}^{b} f(x) dx = -\frac{h^2}{12} (b - a) f^{\prime\prime}(c)
//! ```
//!
//! for some point $c$ where $a ≤ c ≤ b$.
//!
//! A corollary of which is that if $f^{\prime\prime}(x) = 0$ for all $x$ in $\[a,b\]$, i.e. if $f(x)$ is linear,
//! then the trapezoidal rule is exact.
//!
//! The Euler-Maclaurin summation formula also shows that usually n should be chosen large enough
//! so that
//!
//! ```math
//! h = \frac{b-a}{n} < 1
//! ```
//!
//! For example, if $h = 0.1$ then
//! ```math
//! T_{0.1}(f) = \int_{a}^{b} f(x) dx  - 0.00083 \left[f^\prime(b) - f^\prime(a) \right] + 0.00000014 \left[f^{(3)}(b) - f^{(3)}(a) \right] + ...
//! ```
//! and if $h = 0.01$ then
//! ```math
//! T_{0.01}(f) = \int_{a}^{b} f(x) dx  - 0.0000083 \left[ f^\prime(b) - f^\prime(a) \right] + 0.000000000014 \left[f^{(3)}(b) - f^{(3)}(a)\right] + ...
//! ```
//! while if $h=10$ then
//! ```math
//! T_{10}(f) = \int_{a}^{b} f(x) dx  - 8.3333 \left[ f^\prime(b) - f^\prime(a) \right] + 13.89 \left[ f^{(3)}(b) - f^{(3)}(a) \right] + ...
//! ```
//! However, if the function $f(x)$ is linear, then $n$ may be chosen to be $1$.

// extern crate test;

use num::{Float, ToPrimitive, Unsigned};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use super::utils::check_newton_method_args;

/// This function integrates $f(x)$ from $a$ to $a+nh$ using the Simpson's
/// rule by summing from the left end of the interval to the right end.
///  
/// * `f` - Integrand function of a single variable.
/// * `a` - lower limit of the integration interval.
/// * `b` - upper limit of the integration interval.
/// * `n` - number of subintervals.
///
/// # Examples
/// ```
/// use integrator::newton_cotes::trapezoidal::trapezoidal_rule;
///
///
/// fn square(x: f64) -> f64 {
///     x.powi(2)
/// }
///
/// let a = 0.0;
/// let b = 1.0;
///
/// let num_steps: usize = 1_000_000;
///
/// let integral = trapezoidal_rule(square, a, b, num_steps);
/// ```
///
/// # Resources
/// [Methods of numerical Integration (2nd edition), by Philip J. Davis and Philip Rabinowitz.](https://www.cambridge.org/core/journals/mathematical-gazette/article/abs/methods-of-numerical-integration-2nd-edition-by-philip-j-davis-and-philip-rabinowitz-pp-612-3650-1984-isbn-0122063600-academic-press/C331158D0392E1D5CD9B0C6ED4EE5F43)
pub fn trapezoidal_rule<F1: Float + Sync, F2: Float + Send, U: Unsigned + ToPrimitive + Copy>(
    f: fn(F1) -> F2,
    a: F1,
    b: F1,
    n: U,
) -> f64 {
    // checking arguments
    check_newton_method_args(a, b, n);

    // length of each subinterval
    let h: F1 = (b - a) / F1::from(n).expect("failed to convert length of subinterval h");

    // first term of the sum
    let i_0 = f(a).to_f64().unwrap();

    let integral: f64 = (1..(n.to_usize().unwrap()))
        .into_par_iter()
        .map(|i| {
            // subinterval index (as real)
            let i = F1::from(i).expect("failed to convert subinterval index i");
            f(a + i * h).to_f64().unwrap()
        })
        .sum();

    let n: F1 = F1::from(n).expect("failed to convert number of steps n");
    // last term of the sum
    let i_n = f(a + h * n).to_f64().unwrap();

    (0.5 * i_0 + integral + 0.5 * i_n) * h.to_f64().expect("failed to convert subintervql length")
}

#[cfg(test)]
mod tests {
    use std::ops::Div;

    use super::*;
    // use test::Bencher;

    const EPSILON: f64 = 10e-7;
    const NUM_STEPS: usize = 1_000_000;

    #[test]
    fn test_integral_value() {
        fn square(x: f64) -> f64 {
            x.powi(2)
        }

        let a = 0.0;
        let b = 1.0;

        let integral = trapezoidal_rule(square, a, b, NUM_STEPS);

        let analytic_result: f64 = 1.0.div(3.0);

        assert!((integral - analytic_result).abs() < EPSILON);
    }

    #[test]
    fn test_f32_to_f64() {
        // f32 to f64
        fn square(x: f32) -> f64 {
            x.powi(2) as f64
        }

        let a = 0.0;
        let b = 1.0;

        let integral = trapezoidal_rule(square, a, b, NUM_STEPS);

        let analytic_result: f64 = 1.0.div(3.0);

        assert!((integral - analytic_result).abs() < EPSILON);
    }

    #[test]
    fn test_f64_to_f32() {
        // f64 to f32
        fn square(x: f64) -> f32 {
            x.powi(2) as f32
        }

        let a = 0.0;
        let b = 1.0;

        let integral = trapezoidal_rule(square, a, b, NUM_STEPS);

        let analytic_result: f64 = 1.0.div(3.0);

        assert!((integral - analytic_result).abs() < EPSILON);
    }

    #[test]
    fn test_f32_to_f32() {
        // f32 to f32
        fn square(x: f32) -> f32 {
            x.powi(2)
        }

        let a = 0.0;
        let b = 1.0;

        let integral = trapezoidal_rule(square, a, b, NUM_STEPS);

        let analytic_result: f64 = 1.0.div(3.0);

        assert!((integral - analytic_result).abs() < EPSILON);
    }

    // #[bench]
    // fn bench_integral_value(bencher: &mut Bencher) {
    //     fn f1(x: f64) -> f64 {
    //         x.powi(2)
    //     }

    //     let a = 0.0;
    //     let b = 1.0;

    //     bencher.iter(|| {
    //         trapezoidal_rule(f1, a, b, NUM_STEPS);
    //     })
    // }
}
