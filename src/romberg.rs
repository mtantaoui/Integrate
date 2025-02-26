use num::{Float, ToPrimitive, Unsigned};

use rayon::prelude::*;

use std::collections::HashMap;

use std::hash::Hash;

use std::sync::Mutex;

use crate::newton_cotes::trapezoidal_rule;

/// Computes elements of Romberg's matrix recursively given
/// row's and column's index
///
/// * n: Romberg's matrix requested element row index.
/// * m: Romberg's matrix requested element column index.
/// * cache: Storing computed values (shared between threads)
fn romberg<U, F>(
    n: U,
    m: U,
    trapezoids: &[F],
    cache: &Mutex<HashMap<(U, U), F>>, // Shared mutable cache
) -> F
where
    U: Unsigned + ToPrimitive + Send + Copy + Sync + std::hash::Hash + Eq,
    F: Float + Send + Sync,
{
    // Base case
    if m.is_zero() {
        let index = n.to_usize().unwrap();
        return trapezoids[index];
    }

    // Check the cache
    {
        let cache_guard = cache.lock().unwrap();
        if let Some(&value) = cache_guard.get(&(n, m)) {
            return value;
        }
    }

    let one: U = num::one();

    // Compute R[n, m] recursively
    let (r_n_m_minus_1, r_n_1_m_1) = rayon::join(
        || romberg(n, m - one, trapezoids, cache),
        || romberg(n - one, m - one, trapezoids, cache),
    );

    let [coef0, coef1]: [F; 2] = romberg_coefficients(m);
    let result = coef1 * r_n_m_minus_1 - coef0 * r_n_1_m_1;

    // Store in cache
    {
        let mut cache_guard = cache.lock().unwrap();
        cache_guard.insert((n, m), result);
    }

    result
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
/// * `func` - Integrand function of a single variable.
/// * `lower_limit` - lower limit of the integration interval.
/// * `upper_limit` - upper limit of the integration interval.
/// * `n_columns` - number of columns to be used in the Romberg method (columns of the Romberg Matrix).
///
/// This corresponds to a minimum integration subintervals of of length $\dfrac{1}{2^n} * h$
///
///
/// # Resources
/// * [Methods of numerical Integration (2nd edition), by Philip J. Davis and Philip Rabinowitz.](https://www.cambridge.org/core/journals/mathematical-gazette/article/abs/methods-of-numerical-integration-2nd-edition-by-philip-j-davis-and-philip-rabinowitz-pp-612-3650-1984-isbn-0122063600-academic-press/C331158D0392E1D5CD9B0C6ED4EE5F43)
/// * [Romberg's method](https://en.wikipedia.org/wiki/Romberg%27s_method)
pub fn romberg_method<
    Func,
    F1: Float + Sync,
    F2: Float + Sync + Send,
    U: Unsigned + ToPrimitive + Copy + Send + Sync + Hash + Eq,
>(
    func: Func,
    lower_limit: F1,
    upper_limit: F1,
    n_columns: U,
) -> f64
where
    Func: Fn(F1) -> F2 + Sync + Send + Copy,
{
    // first columm of romberg table
    // calculated using trapezoid rule
    let mut trapezoidals: Vec<F2> = Vec::with_capacity(n_columns.to_usize().unwrap());

    // initializing first column of the romberg's matrix using trapezoid rule
    (0..n_columns.to_usize().unwrap())
        .into_par_iter()
        .map(|i| {
            let pow_2 = 2_usize.pow(i.try_into().unwrap()); // 2 ** i
            let trapezoidal = trapezoidal_rule(func, lower_limit, upper_limit, pow_2);
            F2::from(trapezoidal).unwrap()
        })
        .collect_into_vec(&mut trapezoidals);

    // Storing computed values (shared between threads)
    let cache: Mutex<HashMap<(U, U), F2>> = Mutex::new(HashMap::new());

    let integral = romberg(
        n_columns - num::one(),
        n_columns - num::one(),
        trapezoidals.as_slice(),
        &cache,
    );

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
