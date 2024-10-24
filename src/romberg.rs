//! Romberg's method
//!
//! Romberg's method is used to estimate the integral of a function on a closed and bounded interval.
//! Classically, the method consists of successively applying the composite trapezoidal rule, each time
//! halving the length of the subintervals, and using a linear combination of the resulting sequence of
//! estimates to estimate the integral, by successively deleting the low order error terms in
//! the Euler-Maclaurin summation formula.
//! The process terminates when the change of the estimate is within a preassigned tolerance, within
//! a preassigned number of successive estimates.
//!
//! The Euler-Maclaurin summation formula relates the integral of a function $f(x)$ over
//! a closed and bounded interval $\[a,b\]$ , $\int_{a}^{b} f(x) dx$, and the composite trapezoidal rule,
//!
//! ```math
//! T_h (f) = h  \left[ \frac{f(a)}{2} + f(a+h) + ··· + f(b-h) + \frac{f(b)}{2} \right]
//! ```
//!
//! by
//!
//! ```math
//! \begin{split}
//! T_h(f) &= \int_{a}^{b} f(x) dx + \left(\frac{h^2}{12}\right) \left[f^\prime(b) - f^\prime(a)\right] - \left(\frac{h^4}{720}\right) \left[f^{(3)}(b) - f^{(3)}(a) \right] \\
//! &+ ... + K h^{2p-2} \left[ f^{(2p-3)}(b) - f^{(2p-3)(a)} \right] + O(h^{2p})
//! \end{split}
//! ```
//!
//!
//! where $f^\prime$, $f^{(3)}$, and $f^{(2p-3)}$ are the first, third and $(p-3)rd$ derivatives
//! of $f$ and $K$ is a constant.
//!
//! If the subinterval length is halved, then
//!
//! ```math
//! \begin{split}
//! T_{\frac{h}{2}}(f) &= \int_{a}^{b} f(x) dx + \frac{h^2}{4·12} \left[ f^\prime(b) - f^\prime(a) \right] - \left( \frac{h^4}{16·720} \right) \left[ f^{(3)}(b) - f^{(3)}(a) \right] \\
//! &+ ... + K \left( \frac{h}{2} \right)^{2p-2}  \left[ f^{(2p-3)}(b) - f^{(2p-3)(a)} \right] + O(h^{2p})
//! \end{split}
//! ```
//!
//! So that
//!
//! ```math
//! \begin{split}
//! \frac{4 T_{\frac{h}{2}}(f) - T_{h}(f)}{3} &= \int_{a}^{b} f(x) dx + \left( \frac{h^4}{2880} \right) \left[ f^{(3)}(b) - f^{(3)}(a) \right] - \left( \frac{h^6}{96768} \right) \left[ f^{(5)}(b) - f^{(5)}(a) \right] \\
//! &+ ... + K h^{2p-2} \left[ f^{(2p-3)}(b) - f^{(2p-3)(a)} \right] + O(h^{2p})
//! \end{split}
//! ```
//!
//!
//! The $h^2$ term has vanished. This process can be continued, each halving of the subinterval
//! length results in a new composite trapezoidal rule estimate of the integral, which can be
//! combined with previous estimates to yield an estimate, in which the lowest order term
//! involving $h$ vanishes.
//!
//! The easiest way to combine all the estimates from applications
//! of the trapezoidal rule by halving the length of the subintervals is to arrange the estimates
//! in a column, from the coarsest estimate to the finest estimate. To form the second,
//! column take two adjacent values from the first column, subtract the finer estimate
//! from the coarser estimate divide by 3 and add to the finer estimate.
//!
//! The second column is automatically arranged from the coarsest estimate to the finest
//! with one less element. Continue, form the third column by taking two adjacent values
//! from the second column, subtract the finer estimate from the coarser estimate,
//! divide by 15 and add to the finer estimate.
//!
//! The third column is automatically arranged from the coarsest to the finest estimate with one fewer element than in the
//! second column.
//!
//! This process is continued until there is only one element in the last column, this
//! is the estimate of the integral.
//!
//! The numbers which are used the divide the difference of two adjacent elements in the $i^{th}$ column is $4^i - 1$.

use num::{Float, ToPrimitive, Unsigned};

use rayon::prelude::*;

use crate::newton_cotes::trapezoidal::trapezoidal_rule;

/// Computes elements of Romberg's matrix recursively given
/// row's and column's index
///
/// * n: Romberg's matrix requested element row index.
/// * m: Romberg's matrix requested element column index.
fn romberg<U: Unsigned + ToPrimitive + Send + Copy + Sync, F: Float + Send + Sync>(
    n: U,
    m: U,
    trapezoids: &[F],
) -> F {
    if m.is_zero() {
        let index = n.to_usize().unwrap();
        return trapezoids[index];
    }

    let one: U = num::one();

    // R[i,j]: element of Romberg's matrix
    // r_n_m_minus_1: R[n, m-1]
    // r_n_1_m_1: R[n-1, m-1]
    let (r_n_m_minus_1, r_n_1_m_1) = rayon::join(
        || romberg(n, m - one, trapezoids),
        || romberg(n - one, m - one, trapezoids),
    );

    let [coef0, coef1]: [F; 2] = romberg_coefficients(m);
    coef1 * r_n_m_minus_1 - coef0 * r_n_1_m_1
}

/// Approximates the integral of $f(x)$ on $\left[ a, b \right]$ using $T_h(f)$.
///
/// If $T_h(f)$ is the result of applying the trapezoidal rule to approximating
/// the integral of $f(x)$ on $\[a, b\]$ using subintervals of length $h$,   
/// and if $\int_{a}^{b} f(x) dx$ is the integral of $f(x)$ on $\[a,b\]$, then
/// ```math
/// \int_{a}^{b} f(x) dx = \lim_{h \to 0}  T_h(f)
/// ```
/// where the limit is taken as h approaches 0.    
///                          
/// The classical Romberg method applies Richardson Extrapolation to the
/// limit of the sequence $T_h(f), T_{\frac{h}{2}}(f), T_{\frac{h}{4}}(f), ... ,$
/// in which the limit is approached by successively deleting error terms
/// in the Euler-MacLaurin summation formula.  
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
/// # Inputs
/// * `f` - Integrand function of a single variable.
/// * `a` - lower limit of the integration interval.
/// * `b` - lower limit of the integration interval.
/// * `n` - number of columns to be used in the Romberg method (columns of the Romberg Matrix).
///
/// This corresponds to a minimum integration subintervals of of length $\dfrac{1}{2^n} * h$
///
///
/// # Resources
/// * [Methods of numerical Integration (2nd edition), by Philip J. Davis and Philip Rabinowitz.](https://www.cambridge.org/core/journals/mathematical-gazette/article/abs/methods-of-numerical-integration-2nd-edition-by-philip-j-davis-and-philip-rabinowitz-pp-612-3650-1984-isbn-0122063600-academic-press/C331158D0392E1D5CD9B0C6ED4EE5F43)
/// * [Romberg's method](https://en.wikipedia.org/wiki/Romberg%27s_method)
pub fn romberg_method<
    F1: Float + Sync,
    F2: Float + Sync + Send,
    U: Unsigned + ToPrimitive + Copy + Send + Sync,
>(
    f: fn(F1) -> F2,
    a: F1,
    b: F1,
    n: U,
) -> f64 {
    // first columm of romberg table
    // calculated using trapezoid rule
    let mut trapezoidals: Vec<F2> = Vec::with_capacity(n.to_usize().unwrap());

    // initializing first column of the romberg's matrix using trapezoid rule
    (0..n.to_usize().unwrap())
        .into_par_iter()
        .map(|i| {
            let pow_2 = 2_usize.pow(i.try_into().unwrap()); // 2 ** i
            let trapezoidal = trapezoidal_rule(f, a, b, pow_2);
            F2::from(trapezoidal).unwrap()
        })
        .collect_into_vec(&mut trapezoidals);

    let integral = romberg(n - num::one(), n - num::one(), trapezoidals.as_slice());
    integral.to_f64().unwrap()
}

/// Returns coefficients to be used in the Richardson extrapolation for computing
/// Romberg's matrix elements
/// * `m` - order of convergence of Richardson extrapolation.
fn romberg_coefficients<F: Float, U: Unsigned + ToPrimitive>(m: U) -> [F; 2] {
    let m = m.to_i32().unwrap();

    let one = F::from(1.0).unwrap(); // 1

    let _4_m = F::from(4.0.powi(m)).unwrap(); // 4^m
    let _4_m_minus_1 = F::from(4.0.powi(m) - 1.0).unwrap(); // 4^m - 1

    let denominator = one.div(_4_m_minus_1); // 1 / (4^m - 1)

    [
        denominator,        // 1 / (4^m - 1)
        _4_m * denominator, // 4^m / (4^m - 1)
    ]
}

#[cfg(test)]
mod tests {
    use std::ops::Div;

    use super::*;
    // use test::Bencher;

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
        // f32 to f64
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
        // f64 to f32
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
        // f32 to f32
        fn square(x: f32) -> f32 {
            x.powi(2)
        }

        let a = 0.0;
        let b = 1.0;

        let integral = romberg_method(square, a, b, NUM_STEPS);

        let analytic_result: f64 = 1.0.div(3.0);

        assert!((integral - analytic_result).abs() < EPSILON);
    }
}
