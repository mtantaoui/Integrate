use std::{fmt::Debug, iter::Sum, ops::AddAssign};

use num::{bigint::ToBigInt, Float, ToPrimitive, Unsigned};
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};

use crate::utils::{
    chebyshev::{roots_first_kind_chebyshev, roots_second_kind_chebyshev},
    checkers::check_gauss_rule_args,
    hermite::roots_hermite,
    laguerre::roots_laguerre,
    legendre::glpair,
};

/// Approximate the integral of $f(x)$ from $a$ to $b$ using the n point Gauss-Legendre integral approximation formula.
/// * `lower_limit` - lower limit of integration.
/// * `upper_limit` - upper limit of integration.
/// * `n` - number of points to use for Gauss-Legendre integral approximation formula.
///
/// # Examples
///
/// Descriptive of the example here.
///
/// ```
/// use integrate::gauss_quadrature::legendre_rule;
///
///
/// let square = |x: f64| x * x;
///
/// let a = 0.0;
/// let b = 1.0;
///
/// let num_steps: usize = 1_000_000;
///
/// let integral = legendre_rule(square, a, b, num_steps);
/// ```
pub fn legendre_rule<
    Func,
    F1: Float + Sync,
    F2: Float,
    U: Unsigned + ToPrimitive + Copy + PartialOrd + Sync,
>(
    func: Func,
    lower_limit: F1,
    upper_limit: F1,
    n: U,
) -> f64
where
    Func: Fn(F1) -> F2 + Sync,
{
    let two = F1::one() + F1::one();

    let c = (upper_limit - lower_limit) / two;
    let d = (upper_limit + lower_limit) / two;

    let n = n.to_usize().unwrap();

    let integral: f64 = (1..=n)
        .into_par_iter()
        .map(|k| {
            // getting gauss-legendre weights and nodes
            let (_, weight, x) = glpair(n, k);

            // converting node to F1
            let x = F1::from(x).unwrap();

            // interval change formula
            weight * func(c * x + d).to_f64().unwrap() * c.to_f64().unwrap()
        })
        .sum();
    integral
}

/// Approximate the integral of $f(x) e^{-x}$ from 0 to infinity using the $n$
/// point Gauss-Laguerre integral approximation formula.
///
/// * `func` - Integrand function of a single variable.
/// * `n` -  order, number of points used in the rule.  
///
/// # Examples
/// ```
/// use integrate::gauss_quadrature::gauss_laguerre_rule;
///
/// let f = |x: f64| 1.0;
///
/// let n:usize = 100;
///
/// let integral = gauss_laguerre_rule(f, n);
/// ```
pub fn gauss_laguerre_rule<Func, F: Float + Debug + Sync + Send + AddAssign + Sum>(
    func: Func,
    n: usize,
) -> F
where
    Func: Fn(F) -> F + Sync,
{
    check_gauss_rule_args(n);
    let (zeros, weights) = roots_laguerre::<F>(n);

    weights
        .into_par_iter()
        .zip(zeros)
        .map(|(w, x)| w * func(x))
        .sum()
}

/// Approximate the integral of $f(x) e^{-x^2}$ from $-\infty$ to $+\infty$
/// using the $n$ point Gauss-Hermite integral approximation formula.
///
/// The n-th Hermite polynomial is
///
/// $$ H_n(x) = (-1)^n * e^{x^2}* \frac{\partial^{n} e^{-x^2}}{\partial x^n} $$
///            
/// For the n point Gauss-Hermite integral approximation formula the           
/// coefficients are:
///
/// $$A_i = \frac{2^{n+1} * n! * \sqrt{\pi}}{H_{n-1} (x_i)^2} $$
///
/// where $x_i$ is a zero of the n-th Hermite polynomial $H_n(x)$.
///
/// Note that if $x$ is a zero of $H_n(x)$ then $-x$ is also a zero of $H_n(x)$ and the
/// coefficients associated with $x$ and $-x$ are equal.
///
/// # Arguments
///
/// * `func` - Integrand function of a single variable.
/// * `n` -  order, number of points used in the rule.  
///
/// # Examples
/// ```
/// use integrate::gauss_quadrature::gauss_hermite_rule;
///
/// let f = |x: f64| 1.0;
///
/// let n:usize = 100;
///
/// let integral = gauss_hermite_rule(f, n);
/// ```
pub fn gauss_hermite_rule<Func, F: Float + Debug + Sync + Send + AddAssign + Sum + ToBigInt>(
    func: Func,
    n: usize,
) -> F
where
    Func: Fn(F) -> F + Sync,
{
    check_gauss_rule_args(n);

    let (zeros, weights) = roots_hermite::<F>(n);

    weights
        .into_par_iter()
        .zip(zeros)
        .map(|(w, x)| w * func(x))
        .sum()
}

/// Approximate the integral of $\frac{f(x)}{\sqrt{1 - x^2}}$ from -1 to 1
/// using the $n$ point Gauss-Chebyshev first kind integral approximation formula.
///
/// * `func` - Integrand function of a single variable.
/// * `n` -  order, number of points used in the rule.
///
/// # Examples
/// ```
/// use integrate::gauss_quadrature::gauss_first_kind_chebyshev_rule;
///
/// let f = |x: f64| 1.0;
///
/// let n:usize = 100;
///
/// let integral = gauss_first_kind_chebyshev_rule(f, n);
pub fn gauss_first_kind_chebyshev_rule<Func, F: Float + Debug + Sync + Send + AddAssign + Sum>(
    func: Func,
    n: usize,
) -> F
where
    Func: Fn(F) -> F + Sync,
{
    check_gauss_rule_args(n);

    let (zeros, weights) = roots_first_kind_chebyshev::<F>(n);

    weights
        .into_par_iter()
        .zip(zeros)
        .map(|(w, x)| w * func(x))
        .sum()
}

/// Approximate the integral of $f(x) * \sqrt{1 - x^2}$ from -1 to 1
/// using the $n$ point Gauss-Chebyshev second kind integral approximation formula.
///
/// * `func` - Integrand function of a single variable.
/// * `n` -  order, number of points used in the rule.
///
/// # Examples
/// ```
/// use integrate::gauss_quadrature::gauss_second_kind_chebyshev_rule;
///
/// fn f(x: f64) -> f64 {
///     1.0
/// }
///
/// let n:usize = 100;
///
/// let integral = gauss_second_kind_chebyshev_rule(f, n);
pub fn gauss_second_kind_chebyshev_rule<Func, F: Float + Debug + Sync + Send + AddAssign + Sum>(
    f: Func,
    n: usize,
) -> F
where
    Func: Fn(F) -> F + Sync,
{
    check_gauss_rule_args(n);

    let (zeros, weights) = roots_second_kind_chebyshev::<F>(n);

    weights
        .into_par_iter()
        .zip(zeros)
        .map(|(w, x)| w * f(x))
        .sum()
}

#[cfg(test)]
mod legendre_quadrature_tests {
    use std::ops::Div;

    use super::*;

    const EPSILON: f64 = 10e-5;
    const NUM_STEPS: usize = 1_000_000;

    /// Test the numerical integration of ln(x) over the range ]0,1],
    /// exact value of the numerical integration is -1.
    /// Normally, one would not use Gauss-Legendre quadrature for this,
    /// but for the sake of having an example with l > 100, this is included.
    #[test]
    fn test_legendre_rule() {
        fn square(x: f64) -> f64 {
            x.powi(2)
        }

        let a = 0.0;
        let b = 1.0;

        let integral = legendre_rule(square, a, b, NUM_STEPS);

        let analytic_result: f64 = 1.0.div(3.0);

        assert!((integral - analytic_result).abs() < EPSILON);
    }
}

#[cfg(test)]
mod tests {
    use crate::gauss_quadrature::{
        gauss_first_kind_chebyshev_rule, gauss_second_kind_chebyshev_rule,
    };

    const EPSILON: f64 = 10e-6;

    // Test the numerical integration of cos(1000 x) over the range [-1,1]
    // for varying number of Gauss-Chebyshev First Kind quadrature nodes l.
    // exact value of the numerical integration is 0.002 * sin(1000)
    // The fact that only twelve digits of accuracy are obtained is due to the
    // condition number of the summation.
    #[test]
    fn test_chebyshev_first_kind_rule() {
        let exact: f64 = 0.002 * (1000.0_f64).sin();

        fn f(x: f64) -> f64 {
            (1000.0 * x).cos() * (1.0 - x.powi(2)).sqrt()
        }

        println!("Integral Exact Value: {}", exact);

        for l in (540..=700_usize).step_by(20) {
            // Gauss-Legendre rule using glpair function
            let integral: f64 = gauss_first_kind_chebyshev_rule(f, l);

            println!(
                "number of nodes: {} \t Gauss-Chebyshev First Kind Integral: {}",
                l, integral
            );

            assert!(integral - exact < EPSILON);
        }
    }

    // Test the numerical integration of cos(1000 x) over the range [-1,1]
    // for varying number of Gauss-Chebyshev Second Kind quadrature nodes l.
    // exact value of the numerical integration is 0.002 * sin(1000)
    // The fact that only twelve digits of accuracy are obtained is due to the
    // condition number of the summation.
    #[test]
    fn test_chebyshev_second_kind_rule() {
        let exact: f64 = 0.002 * (1000.0_f64).sin();

        fn f(x: f64) -> f64 {
            (1000.0 * x).cos() / (1.0 - x.powi(2)).sqrt()
        }

        println!("Integral Exact Value: {}", exact);

        for l in (540..=700_usize).step_by(20) {
            // Gauss-Legendre rule using glpair function
            let integral: f64 = gauss_second_kind_chebyshev_rule(f, l);

            println!(
                "number of nodes: {} \t Gauss-Chebyshev Second Kind Integral: {}",
                l, integral
            );

            assert!((integral - exact).abs() < EPSILON);
        }
    }
}
