//! Rectangle Rule
//!
//! The rectangle rule approximates the integral of a function $f(x)$ on the
//! closed and bounded interval $[a, a+h]$ of length $h > 0$ by the (signed) area
//! of the rectangle with length h and height the value of the function $f(x)$
//! evaluated at the midpoint of the interval, $f(a+h/2)$.
//!
//! The composite rectangle rule is used to approximate the integral of a function
//! $f(x)$ over a closed and bounded interval $[a, b]$ where $a < b$, by decomposing
//! the interval $[a, b]$ into $n > 1$ subintervals of equal length $h = \frac{b-a}{n}$
//! and adding the results of applying the rectangle rule to each subinterval.
//!
//! By abuse of language both the composite rectangle rule and the rectangle rule sometimes
//! are referred to simply as the rectangle rule.
//!
//! Let $\int_{a}^{b} f(x) dx$ be the integral of f(x) over the closed and bounded interval $\[a ,b \]$,
//! and let $R_h(f)$ be the result of applying the rectangle rule with n subintervals of length h, i.e.
//!
//! ```math
//! R_h(f)=h \left[ f(a+\frac{h}{2}) + f(a+\frac{3h}{2}) + ··· + f(b-\frac{h}{2}) \right]
//! ```
//!
//! An immediate consequence of the Euler-Maclaurin summation formula yields the following equation
//! relating $\int_{a}^{b} f(x) dx$ and $R_h(f)$:
//!
//! ```math
//! \begin{split}
//! R_h(f) & = \int_{a}^{b} f(x) dx - \frac{h^2}{24} \left[ f^\prime (b) - f^\prime (a) \right] +  \frac{7h^4}{5760} \left[ f^{(3)}(b) - f^{(3)}(a) \right] + \\
//! & ··· + K h^{2p-2} \left[ f^{(2p-3)}(b) - f^{(2p-3)}(a) \right] + O(h^{2p})  
//! \end{split}
//! ```
//!
//! where $f'$, $f^{(3)}$, and $f^{(2p-3)}$ are the first, third and $(2p-3)^{rd}$ derivatives of $f$ and $K$ is a constant.
//!
//!
//! The last term, $O(h^{2p})$ is important. Given an infinitely differentiable function
//! in which the first $2p-3$ derivatives vanish at both endpoints of the interval of integration,
//! it is not true that $R_h(f) = \int_{a}^{b} f(x) dx$, but rather what the theorem says is that
//!
//! ```math
//! \lim_{h \to 0} \left| \dfrac{R_h(f) - \int_{a}^{b} f(x)dx}{h^{2p}} \right| < M
//! ```
//!
//! where $M>0$.
//!
//! If $f$ is at least twice differentiable on the interval $\[a,b\]$, then applying the mean-value
//! theorem to
//!
//! ```math
//! \begin{split}
//! R_h(f) - \int_{a}^{b} f(x) dx & = -\frac{h^2}{24} \left[ f^\prime (b) - f^\prime (a) \right] +  \frac{7h^4}{5760} \left[ f^{(3)}(b) - f^{(3)}(a) \right] \\
//! & +  ··· + K h^{2p-2} \left[ f^{(2p-3)}(b) - f^{(2p-3)}(a) \right] + O(h^{2p})
//! \end{split}
//! ```
//!
//! yields the standard truncation error expression
//!
//! ```math
//! R_h(f) - \int_{a}^{b} f(x) dx = -\frac{h^2}{24} (b - a) f^{\prime\prime}(c)
//! ```
//!
//! for some point $c$ where $a ≤ c ≤ b$.
//!
//! A corollary of which is that if $f''(x) = 0$ for all $x$ in $\[a,b\]$,
//! i.e. if $f(x)$ is linear, then the rectangle rule is exact.
//!
//! The Euler-Maclaurin summation formula also shows that usually $n$ should be chosen large enough
//! so that $h = (b - a) / n < 1$. For example, if h = 0.1 then
//!
//! ```math
//! \begin{split}
//! R_{0.1}(f) &= \int_{a}^{b} f(x) dx  - 0.00042 \left[ f'(b) - f'(a) \right]  \\
//! & + 0.00000012 \left[f^{(3)}(b) - f^{(3)}(a) \right] + ...
//! \end{split}
//! ```
//! and if $h = 0.01$ then
//!
//! ```math
//! \begin{split}
//! R_{0.01}(f) &= \int_{a}^{b} f(x) dx  - 0.0000042 \left[ f^\prime(b) - f^\prime(a) \right] \\
//! &+ 0.000000000012 \left[ f^{(3)}(b) - f^{(3)}(a) \right] + ...   
//! \end{split}
//! ```
//!
//! while if $h=10$ then
//! ```math
//! \begin{split}
//! R_{10}(f) &= \int_{a}^{b} f(x) dx  - 4.1667 \left[ f^\prime(b) - f^\prime(a)\right] \\
//! &+ 12.15 \left[ f^{(3)}(b) - f^{(3)}(a) \right] + ...
//! \end{split}
//! ```
//! However, if the function $f(x)$ is linear, then $n$ may be chosen to be $1$.

// extern crate test;

use num::{Float, ToPrimitive, Unsigned};
use rayon::prelude::*;

use super::utils::check_newton_method_args;

/// This function integrates $f(x)$ from $a$ to $a+nh$ using the rectangle
/// rule by summing from the left end of the interval to the right end.
///
/// * `f` - Integrand function of a single variable.
/// * `a` - lower limit of the integration interval.
/// * `b` - upper limit of the integration interval.
/// * `n` - number of subintervals.
///
/// # Examples
/// ```
/// use integrate::newton_cotes::rectangle::rectangle_rule;
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
/// let integral = rectangle_rule(square, a, b, num_steps);
/// ```
///
/// # Resources
/// [Methods of numerical Integration (2nd edition), by Philip J. Davis and Philip Rabinowitz.](https://www.cambridge.org/core/journals/mathematical-gazette/article/abs/methods-of-numerical-integration-2nd-edition-by-philip-j-davis-and-philip-rabinowitz-pp-612-3650-1984-isbn-0122063600-academic-press/C331158D0392E1D5CD9B0C6ED4EE5F43)
pub fn rectangle_rule<F1: Float + Sync, F2: Float, U: Unsigned + ToPrimitive + Copy>(
    f: fn(F1) -> F2,
    a: F1,
    b: F1,
    n: U,
) -> f64 {
    // checking arguments
    check_newton_method_args(a, b, n);

    // length of each subinterval
    let h: F1 = (b - a) / F1::from(n).expect("failed to convert length of subinterval h");

    let integral: f64 = (0..(n.to_usize().unwrap()))
        .into_par_iter()
        .map(|i| {
            // subinterval index (as real)
            let i = F1::from(i).expect("failed to convert subinterval index i");

            // subinterval midpoint
            let x = a + i * h + (h / F1::from(2).expect("failed to compute subinterval midpoint"));

            // converting f(x) to primitive type f64
            f(x).to_f64().expect("failed to convert f(x) to f64")
        })
        .sum();
    integral * h.to_f64().unwrap()
}

#[cfg(test)]
mod tests {
    use std::ops::Div;

    use super::*;
    // use test::Bencher;

    const EPSILON: f64 = 10e-5;
    const NUM_STEPS: usize = 1000;

    #[test]
    fn test_integral_value() {
        fn square(x: f64) -> f64 {
            x.powi(2)
        }

        let a = 0.0;
        let b = 1.0;

        let integral = rectangle_rule(square, a, b, NUM_STEPS);

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

        let integral = rectangle_rule(square, a, b, NUM_STEPS);

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

        let integral = rectangle_rule(square, a, b, NUM_STEPS);

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

        let integral = rectangle_rule(square, a, b, NUM_STEPS);

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
    //         rectangle_rule(f1, a, b, NUM_STEPS);
    //     })
    // }
}
