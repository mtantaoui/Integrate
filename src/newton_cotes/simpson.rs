//! Simpson's Rule
//!
//! Simpson's rule approximates the integral of a function $f(x)$ on the closed and bounded interval
//! $\[a, a+h\]$ of length $h > 0$ by the integral on $\[a, a+h\]$ of the quadratic passing through the
//! points $\left(a, f(a)\right)$, $\left(a+\dfrac{h}{2}, f(a+\dfrac{h}{2})\right)$ and $\left(a+h, f(a+h)\right)$.
//!
//! The composite Simpson's rule is used to approximate the integral of a function $f(x)$ over a closed and bounded interval $\[a, b\]$ where $a < b$,
//! by decomposing the interval $\[a, b\]$ into $n > 1$ subintervals of equal length $h = \dfrac{b-a}{n}$,
//! then adding the results of applying the Simpson's rule to each subinterval. By abuse of
//! language both the composite Simpson's rule and Simpson's rule sometimes are referred to simply as
//! Simpson's rule. Let $\int_{a}^{b} f(x) dx$ be the integral of $f(x)$ over the closed and bounded interval
//! $\[a,b\]$, and let $S_h(f)$ be the result of applying the Simpson's rule with $n$ subintervals of length h, i.e.
//!
//! ```math
//! \begin{split}
//! S_h(f) &= \frac{h}{6} \left[ f(a) + 4f(a+\frac{h}{2}) + 2f(a+h) \right. \\
//! & + \left. ··· + 2f(b-h) + 4f(b-\frac{h}{2}) + f(b) \right]
//! \end{split}
//! ```
//!
//! An immediate consequence of the Euler-Maclaurin summation formula relates $\int_{a}^{b} f(x) dx$ and $S_h(f)$
//! ```math
//! \begin{split}
//! S_h(f) & = \int_{a}^{b} f(x) dx + \frac{h^4}{2880} \left[ f^{(3)}(b) - f^{(3)}(a) \right] - \frac{h^6}{96768} \left[ f^{(5)}(b) - f^{(5)}(a) \right] \\
//! & + ··· + K h^{2p-2} \left[ f^{(2p-3)}(b) - f^{(2p-3)}(a) \right] + O(h^{2p})
//! \end{split}
//! ```
//! where $f^{(3)}$, $f^{(5)}$, and $f^{(2p-3)}$ are the third, fifth and $(p-3)^{th}$ derivatives of $f$ and $K$ is a constant.
//!
//! The last term, $O(h^{2p})$ is important. Given an infinitely differentiable function in which the first
//! $2p-3$ derivatives vanish at both endpoints of the interval of integration, it is not true that
//! ```math
//! S_h(f) = \int_{a}^{b}f( x ) dx
//! ```
//! but rather what the theorem says is that
//! ```math
//! \lim_{h \to 0} \left| \frac{S_h(f) - \int_{a}^{b} f(x)dx}{h^{2p}} \right| < M
//! ```
//! where M > 0.
//!
//! If $f$ is at least four times differentiable on the interval $\[a,b\]$, then applying the mean-value theorem to
//! ```math
//! \begin{split}
//! S_h(f) - \int_{a}^{b}f( x ) dx &= \frac{h^4}{2880} \left[ f^{(3)}(b) - f^{(3)}(a) \right] - \frac{h^6}{96768} \left[ f^{(5)}(b) - f^{(5)}(a) \right] \\
//! &+ ··· + K h^{(2p - 2)} \left[f^{(2p-3)}(b) - f^{(2p-3)}(a) \right] + O(h^{2p})
//! \end{split}
//! ```
//! yields the standard truncation error expression
//!
//! ```math
//! S_h(f) - \int_{a}^{b} f( x ) dx = \frac{h^4}{2880} (b-a) f^{(4)}(c)
//! ```
//!  for some point $c$ where $a ≤ c ≤ b$.
//!
//!
//! A corollary of which is that if $f^{(4)}(x) = 0$ for all $x$ in $\[a,b\]$, i.e. if $f(x)$ is a cubic,
//! then Simpson's rule is exact.
//! The Euler-Maclaurin summation formula also shows that usually $n$ should be chosen large enough so that $h = \frac{b-a}{n} < 1$.
//!
//! For example, if $h = 0.1$ then
//!
//! ```math
//! \begin{split}
//! S_{0.1}(f) & = \int_{a}^{b} f(x) dx + 3.5 · 10^{-8} \left[ f^{3}(b) - f^{3}(a) \right]\\
//!  & - 1.033 · 10^{-11} \left[ f^{(5)}(b) - f^{(5)}(a) \right] + ···
//! \end{split}
//! ```
//!
//!
//! and if $h = 0.01$ then
//!
//! ```math
//! \begin{split}
//! S_{0.01}(f) & = \int_{a}^{b} f(x) dx + 3.5 · 10^{-12} \left[ f^{3}(b) - f^{3}(a) \right] \\
//! & - 1.033 · 10^{-17} \left[ f^{(5)}(b) - f^{(5)}(a) \right] + ···
//! \end{split}
//! ```
//!
//! while if $h = 10$ then
//!
//! ```math
//! \begin{split}
//! S_{0.01}(f) & = \int_{a}^{b} f(x) dx + 3.47  \left[ f^{3}(b) - f^{3}(a) \right] \\
//! & - 10.33  \left[ f^{(5)}(b) - f^{(5)}(a) \right] + ···
//! \end{split}
//! ```
//!
//! However, if the function $f(x)$ is a cubic, then $n$ may be chosen to be $1$.

use std::ops::Div;

use num::{Float, ToPrimitive, Unsigned};

use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};

use super::utils::check_newton_method_args;

/// This function integrates $f(x)$ from $a$ to $a+nh$ using the trapezoidal
/// rule by summing from the left end of the interval to the right end.
///
/// * `f` - Integrand function of a single variable.
/// * `a` - lower limit of the integration interval.
/// * `b` - upper limit of the integration interval.
/// * `n` - number of subintervals.
///
/// # Examples
/// ```
/// use integrate::newton_cotes::simpson::simpson_rule;
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
/// let integral = simpson_rule(square, a, b, num_steps);
/// ```
///
/// # Resources
/// [Methods of numerical Integration (2nd edition), by Philip J. Davis and Philip Rabinowitz.](https://www.cambridge.org/core/journals/mathematical-gazette/article/abs/methods-of-numerical-integration-2nd-edition-by-philip-j-davis-and-philip-rabinowitz-pp-612-3650-1984-isbn-0122063600-academic-press/C331158D0392E1D5CD9B0C6ED4EE5F43)
pub fn simpson_rule<F1: Float + Sync, F2: Float, U: Unsigned + ToPrimitive + Copy>(
    f: fn(F1) -> F2,
    a: F1,
    b: F1,
    n: U,
) -> f64 {
    // checking arguments
    check_newton_method_args(a, b, n);

    // length of each subinterval
    let h: F1 = (b - a) / F1::from(n).expect("failed to convert length of subinterval h");

    // half the length of each subinterval h/2
    let h_over_2 = h / F1::from(2).unwrap();

    // first term of the sum
    let i_0 = f(a).to_f64().unwrap() + 4.0 * f(a + h_over_2).to_f64().unwrap();

    let integral: f64 = (2..(2 * n.to_usize().unwrap()))
        .into_par_iter()
        .step_by(2)
        .map(|i| {
            // subinterval index (as real)
            let i_plus_1 = F1::from(i + 1).expect("failed to convert subinterval index (i+1)");
            let i = F1::from(i).expect("failed to convert subinterval index i");

            2.0 * f(a + i * h_over_2).to_f64().unwrap()
                + 4.0 * f(a + i_plus_1 * h_over_2).to_f64().unwrap()
        })
        .sum();

    let n = F1::from(n).expect("failed to convert n");
    let i_n = f(a + n * h_over_2).to_f64().unwrap();

    (i_0 + integral + i_n) * h.to_f64().unwrap() * 1.0.div(6.0)
}

#[cfg(test)]
mod tests {

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

        let integral = simpson_rule(square, a, b, NUM_STEPS);

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

        let integral = simpson_rule(square, a, b, NUM_STEPS);

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

        let integral = simpson_rule(square, a, b, NUM_STEPS);

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

        let integral = simpson_rule(square, a, b, NUM_STEPS);

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
    //         simpson_rule(f1, a, b, NUM_STEPS);
    //     })
    // }
}
