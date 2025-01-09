//! Newton's 3/8 Rule
//!
//! Newton's 3/8 rule approximates the integral of a function $f(x)$ on the closed and bounded
//! interval $\[a, a+h\]$ of length $h > 0$ by the integral on $\[a, a+h\]$ of the cubic passing
//! through the points $(a, f(a))$, $(a+\frac{h}{3}, f(a+\frac{h}{3}))$, $(a+\frac{2h}{3}, f(a+\frac{2h}{3}))$ and $(a+h, f(a+h))$.
//!
//! The composite Newton's 3/8 rule is used to approximate the integral of a function $f(x)$
//! over a closed and bounded interval $\[a, b\]$ where $a < b$, by decomposing the interval
//! $\[a, b\]$ into $n > 1$ subintervals of equal length $h = \frac{b-a}{n}$, then adding the
//! results of applying the Newton's 3/8 rule to each subinterval.
//!
//! By abuse of language both the composite Newton's 3/8 rule and Newton's 3/8
//! rule are referred to simply as Newton's 3/8 rule. Let $\int_{a}^{b} f(x)dx$ be the
//! integral of $f(x)$ over the closed and bounded interval $\[a, b\]$, and let $N_h(f)$
//! be the result of applying the Newton's 3/8 rule with $n$ subintervals of length $h$, i.e.
//!
//! ```math
//! \begin{split}
//! N_h(f) &=  \frac{h}{8}  \left[ f(a) + 3 f\left(a+ \frac{h}{3} \right) + 3f\left(a+ \frac{2h}{3} \right) + 2 f(a + h) \right.\\
//! & \left. + ··· + 2f(b-h) + 3f \left( b-\frac{2h}{3} \right) + 3f \left( b - \frac{h}{3} \right) + f(b) \right]
//! \end{split}
//! ```
//!
//! An immediate consequence of the Euler-Maclaurin summation formula relates $\int_{a}^{b} f(x)dx$ and $N_h(f)$
//!
//! ```math
//! \begin{split}
//! N_h(f) &= \int_{a}^{b} f(x)dx + \frac{h^4}{6480} \left[ f^{3}(b) - f^{3}(a) \right] - \frac{h^6}{244944} \left[ f^{(5)}(b) - f^{(5)}(a) \right] \\
//! & + ··· + K h^{2p-2} \left[ f^{(2p - 3)}(b) - f^{(2p - 3)}(a) \right] + O(h^{2p})
//! \end{split}
//! ```
//!
//! where $f^{(3)}$, $f^{(5)}$, and $f^{(2p-3)}$ are the third, fifth and $(p-3)rd$ derivatives
//! of $f$ and $K$ is a constant.
//!
//! The last term, $O(h^{2p})$ is important. Given an infinitely differentiable function in which the
//! first $2p-3$ derivatives vanish at both endpoints of the interval of integration, it is not true that
//!
//! ```math
//! \begin{split}
//! N_h(f) = \int_{a}^{b} f(x) dx
//! \end{split}
//! ```
//!
//! but rather what the theorem says is that
//!
//! ```math
//! \begin{split}
//! \lim_{h \to 0} \mid \frac{N_h(f) - \int_{a}^{b} f(x)dx}{h^{2p}} \mid < M
//! \end{split}
//! ```
//!
//! where $M > 0$.
//!
//! If $f$ is at least four times differentiable on the interval $\[a,b\]$, then applying the
//! mean-value theorem to
//!
//! ```math
//! \begin{split}
//! N_h(f) - \int_{a}^{b} f(x)dx &= \frac{h^4}{6480} \left[ f^{(3)}(b) - f^{(3)}(a) \right] - \frac{h^6}{244944} \left[ f^{(5)}(b) - f^{(5)}(a) \right]  \\
//! &  + ··· + K h^{2p - 2} \left[ f^{(2p - 3)}(b) - f^{(2p - 3)}(a) \right] + O(h^{2p})
//! \end{split}
//! ```
//!
//! yields the standard truncation error expression
//!
//! ```math
//! N_h(f) - \int_{a}^{b} f(x)dx = \frac{h^4}{6480}  (b-a) f^{(4)}(c)
//! ```  
//!
//! for some point $c$ where $a ≤ c ≤ b$.
//!
//! A corollary of which is that if $f^{(4)}(x) = 0$ for all $x$ in $\[a,b\]$, i.e. if $f(x)$ is a cubic,
//! then Newton's 3/8 rule is exact.
//!
//! The Euler-Maclaurin summation formula also shows that usually $n$ should be chosen large enough so that $h = \dfrac{b-a}{n}< 1$.
//!
//! For example, if $h = 0.1$ then
//!
//! ```math
//! \begin{split}
//! N_{0.1}(f) &= \int_{a}^{b} f(x)dx  + 1.5·10^{-8} \left[ f^{(3)}(b) - f^{(3)}(a) \right] \\
//!  &- 4.1 · 10^{-12} \left[ f^{(5)}(b) - f^{(5)}(a) \right] + ···
//! \end{split}
//! ```
//!
//! and if $h = 0.01$ then
//!
//! ```math
//! \begin{split}
//! N_{0.01}(f) &= \int_{a}^{b} f(x)dx + 1.5·10^{-12} [ f^{3}(b) - f^{3}(a) ] \\
//! &- 4.1·10^{-18} [ f^{(5)}(b) - f^{(5)}(a) ] + ···
//! \end{split}
//! ```
//! while if $h = 10$ then
//!
//! ```math
//! \begin{split}
//! N_{10}(f) &= \int_{a}^{b} f(x)dx  + 1.54 \left[ f^{(3)}(b) - f^{(3)}(a) \right] \\
//! &- 4.08 \left[ f^{(5)}(b) - f^{(5)}(a) \right] + ···
//! \end{split}
//! ```
//!
//! However, if the function $f(x)$ is a cubic, then $n$ may be chosen to be 1.

// extern crate test;

use std::ops::Div;

use num::{Float, ToPrimitive, Unsigned};

use super::utils::check_newton_method_args;

use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};

/// This function integrates $f(x)$ from $a$ to $a+nh$ using the Newton's 3/8
/// rule by summing from the left end of the interval to the right end.
///
/// * `func` - Integrand function of a single variable.
/// * `lower_limit` - lower limit of the integration interval.
/// * `upper_limit` - upper limit of the integration interval.
/// * `n_intervals` - number of subintervals.
///
/// # Examples
/// ```
/// use integrate::newton_cotes::newton::newton_rule;
///
///
/// let square = |x: f64| x * x;
///
/// let a = 0.0;
/// let b = 1.0;
///
/// let num_steps: usize = 1_000_000;
///
/// let integral = newton_rule(square, a, b, num_steps);
/// ```
///
/// # Resources
/// [Methods of numerical Integration (2nd edition), by Philip J. Davis and Philip Rabinowitz.](https://www.cambridge.org/core/journals/mathematical-gazette/article/abs/methods-of-numerical-integration-2nd-edition-by-philip-j-davis-and-philip-rabinowitz-pp-612-3650-1984-isbn-0122063600-academic-press/C331158D0392E1D5CD9B0C6ED4EE5F43)
pub fn newton_rule<Func, F1: Float + Sync, F2: Float, U: Unsigned + ToPrimitive + Copy>(
    func: Func,
    lower_limit: F1,
    upper_limit: F1,
    n_intervals: U,
) -> f64
where
    Func: Fn(F1) -> F2 + Sync,
{
    // checking arguments
    check_newton_method_args(lower_limit, upper_limit, n_intervals);

    // length of each subinterval
    let h: F1 = (upper_limit - lower_limit)
        / F1::from(n_intervals).expect("failed to convert length of subinterval h");

    // half the length of each subinterval h/3
    let h_over_3 = h / F1::from(3).unwrap();

    // half the length of each subinterval (h * 2)/3
    let h2_over_3 = F1::from(2).unwrap() * (h / F1::from(3).unwrap());

    // first term of the sum
    let i_0 = func(lower_limit).to_f64().unwrap()
        + 3.0
            * (func(lower_limit + h_over_3).to_f64().unwrap()
                + func(lower_limit + h2_over_3).to_f64().unwrap());

    let integral: f64 = (3..(3 * n_intervals.to_usize().unwrap()))
        .into_par_iter()
        .step_by(3)
        .map(|i| {
            // subinterval index (as real)
            let i_plus_1 = F1::from(i + 1).expect("failed to convert subinterval index (i+1)");
            let i_plus_2 = F1::from(i + 2).expect("failed to convert subinterval index (i+2)");
            let i = F1::from(i).expect("failed to convert subinterval index i");

            2.0 * func(lower_limit + i * h_over_3).to_f64().unwrap()
                + 3.0
                    * (func(lower_limit + i_plus_1 * h_over_3).to_f64().unwrap()
                        + func(lower_limit + i_plus_2 * h_over_3).to_f64().unwrap())
        })
        .sum();

    let n = F1::from(n_intervals).expect("failed to convert n");
    let i_n = func(lower_limit + n * h_over_3).to_f64().unwrap();

    (i_0 + integral + i_n) * h.to_f64().unwrap() * 1.0.div(8.0)
}

#[cfg(test)]
mod tests {

    use super::*;
    // use test::Bencher;

    const EPSILON: f64 = 10e-4;
    const NUM_STEPS: usize = 1_000;

    #[test]
    fn test_integral_value() {
        fn square(x: f64) -> f64 {
            x.powi(2)
        }

        let a = 0.0;
        let b = 1.0;

        let integral = newton_rule(square, a, b, NUM_STEPS);

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

        let integral = newton_rule(square, a, b, NUM_STEPS);

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

        let integral = newton_rule(square, a, b, NUM_STEPS);

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

        let integral = newton_rule(square, a, b, NUM_STEPS);

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
    //         newton_rule(f1, a, b, NUM_STEPS);
    //     })
    // }
}
