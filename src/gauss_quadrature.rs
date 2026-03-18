//! Gaussian quadrature rules for numerical integration.
//!
//! Gaussian quadrature approximates $\int f(x)\\,w(x)\\,dx$ by evaluating $f$ at $n$
//! carefully chosen nodes $x_i$ and forming a weighted sum $\sum_{i=1}^n w_i f(x_i)$.
//! The nodes are zeros of the corresponding orthogonal polynomial; the weights are derived
//! from the Golub-Welsch algorithm (eigenvalues of the Jacobi tridiagonal matrix).

use std::{fmt::Debug, iter::Sum, ops::AddAssign};

use num::{bigint::ToBigInt, Float, ToPrimitive, Unsigned};

use crate::utils::{
    chebyshev::{roots_first_kind_chebyshev, roots_second_kind_chebyshev},
    checkers::check_gauss_rule_args,
    hermite::roots_hermite,
    laguerre::roots_laguerre,
    legendre::glpair,
};

/// Approximate $\int_a^b f(x)\\,dx$ using $n$-point Gauss-Legendre quadrature.
///
/// The integral is mapped to $[-1, 1]$ via the linear change of variables
/// $x = \frac{b-a}{2}\\,t + \frac{b+a}{2}$, giving
///
/// $$\int_a^b f(x)\\,dx = \frac{b-a}{2}\int_{-1}^{1} f\left(\tfrac{b-a}{2}t+\tfrac{b+a}{2}\right)dt \approx \frac{b-a}{2}\sum_{i=1}^{n} w_i\\, g(t_i)$$
///
/// where $t_i$ and $w_i$ are the Gauss-Legendre nodes and weights on $[-1, 1]$.
///
/// # Arguments
///
/// * `func` - Integrand $f$, a function of a single variable.
/// * `lower_limit` - Lower limit of integration $a$.
/// * `upper_limit` - Upper limit of integration $b$.
/// * `n` - Number of quadrature points (polynomial order).
///
/// # Examples
///
/// ```
/// use integrate::gauss_quadrature::legendre_rule;
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
///
/// # Resources
///
/// - Davis, P. J. & Rabinowitz, P., *Methods of Numerical Integration*, 2nd ed., Academic Press (1984).
/// - Press et al., *Numerical Recipes in C*, Chapter 4.
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

/// Approximate $\int_0^\infty f(x)\\,e^{-x}\\,dx$ using $n$-point Gauss-Laguerre quadrature.
///
/// $$\int_0^\infty f(x)\\,e^{-x}\\,dx \approx \sum_{i=1}^{n} w_i\\, f(x_i)$$
///
/// where $x_i$ are the zeros of the $n$-th Laguerre polynomial $L_n(x)$ and the
/// weights are
///
/// $$w_i = \frac{x_i}{(n+1)^2 \left[L_{n+1}(x_i)\right]^2}$$
///
/// # Arguments
///
/// * `func` - Integrand $f$, a function of a single variable. The weight $e^{-x}$ is
///   absorbed into the quadrature rule; pass $f$ without it.
/// * `n` - Number of quadrature points (polynomial order).
///
/// # Panics
///
/// Panics if `n == 0` (via internal argument check).
///
/// # Examples
///
/// ```
/// use integrate::gauss_quadrature::gauss_laguerre_rule;
///
/// let f = |x: f64| 1.0;
///
/// let n: usize = 100;
///
/// let integral = gauss_laguerre_rule(f, n);
/// ```
///
/// # Resources
///
/// - [Wikipedia: Gauss-Laguerre quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Laguerre_quadrature)
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
        .into_iter()
        .zip(zeros)
        .map(|(w, x)| w * func(x))
        .sum()
}

/// Approximate $\int_{-\infty}^{+\infty} f(x)\\,e^{-x^2}\\,dx$ using $n$-point Gauss-Hermite quadrature.
///
/// The $n$-th Hermite polynomial is defined as
///
/// $$H_n(x) = (-1)^n\\, e^{x^2}\\, \frac{\partial^n e^{-x^2}}{\partial x^n}$$
///
/// The quadrature weights evaluated at node $x_i$ (a zero of $H_n$) are
///
/// $$A_i = \frac{2^{n+1}\\, n!\\, \sqrt{\pi}}{\left[H_{n-1}(x_i)\right]^2}$$
///
/// giving the approximation
///
/// $$\int_{-\infty}^{+\infty} f(x)\\,e^{-x^2}\\,dx \approx \sum_{i=1}^{n} A_i\\, f(x_i)$$
///
/// Because $H_n$ is either even or odd, if $x$ is a zero then $-x$ is also a zero
/// and both share the same weight.
///
/// > **Note:** For large $n$ some weights may underflow to zero. Check `stderr` for
/// > any warnings emitted by the node/weight computation.
///
/// # Arguments
///
/// * `func` - Integrand $f$, a function of a single variable. The Gaussian weight
///   $e^{-x^2}$ is absorbed into the rule; pass $f$ without it.
/// * `n` - Number of quadrature points (polynomial order).
///
/// # Panics
///
/// Panics if `n == 0` (via internal argument check).
///
/// # Examples
///
/// ```
/// use integrate::gauss_quadrature::gauss_hermite_rule;
///
/// let f = |x: f64| 1.0;
///
/// let n: usize = 100;
///
/// let integral = gauss_hermite_rule(f, n);
/// ```
///
/// # Resources
///
/// - [Wikipedia: Gauss-Hermite quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Hermite_quadrature)
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
        .into_iter()
        .zip(zeros)
        .map(|(w, x)| w * func(x))
        .sum()
}

/// Approximate $\int_{-1}^{1} \frac{f(x)}{\sqrt{1-x^2}}\\,dx$ using $n$-point
/// Gauss-Chebyshev quadrature of the first kind.
///
/// The nodes and weights are
///
/// $$x_i = \cos\left(\frac{(2i-1)\pi}{2n}\right), \qquad w_i = \frac{\pi}{n}$$
///
/// giving the approximation
///
/// $$\int_{-1}^{1} \frac{f(x)}{\sqrt{1-x^2}}\\,dx \approx \frac{\pi}{n}\sum_{i=1}^{n} f(x_i)$$
///
/// All $n$ weights are equal to $\pi/n$.
///
/// # Arguments
///
/// * `func` - Integrand $f$, a function of a single variable. The weight
///   $1/\sqrt{1-x^2}$ is absorbed into the rule; pass $f$ without it.
/// * `n` - Number of quadrature points (polynomial order).
///
/// # Panics
///
/// Panics if `n == 0` (via internal argument check).
///
/// # Examples
///
/// ```
/// use integrate::gauss_quadrature::gauss_first_kind_chebyshev_rule;
///
/// let f = |x: f64| 1.0;
///
/// let n: usize = 100;
///
/// let integral = gauss_first_kind_chebyshev_rule(f, n);
/// ```
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
        .into_iter()
        .zip(zeros)
        .map(|(w, x)| w * func(x))
        .sum()
}

/// Approximate $\int_{-1}^{1} f(x)\sqrt{1-x^2}\\,dx$ using $n$-point
/// Gauss-Chebyshev quadrature of the second kind.
///
/// The nodes and weights are
///
/// $$x_i = \cos\left(\frac{i\pi}{n+1}\right), \qquad w_i = \frac{\pi}{n+1}\sin^2\left(\frac{i\pi}{n+1}\right)$$
///
/// giving the approximation
///
/// $$\int_{-1}^{1} f(x)\sqrt{1-x^2}\\,dx \approx \sum_{i=1}^{n} w_i\\, f(x_i)$$
///
/// # Arguments
///
/// * `func` - Integrand $f$, a function of a single variable. The weight
///   $\sqrt{1-x^2}$ is absorbed into the rule; pass $f$ without it.
/// * `n` - Number of quadrature points (polynomial order).
///
/// # Panics
///
/// Panics if `n == 0` (via internal argument check).
///
/// # Examples
///
/// ```
/// use integrate::gauss_quadrature::gauss_second_kind_chebyshev_rule;
///
/// fn f(x: f64) -> f64 {
///     1.0
/// }
///
/// let n: usize = 100;
///
/// let integral = gauss_second_kind_chebyshev_rule(f, n);
/// ```
pub fn gauss_second_kind_chebyshev_rule<Func, F: Float + Debug + Sync + Send + AddAssign + Sum>(
    f: Func,
    n: usize,
) -> F
where
    Func: Fn(F) -> F + Sync,
{
    check_gauss_rule_args(n);

    let (zeros, weights) = roots_second_kind_chebyshev::<F>(n);

    weights.into_iter().zip(zeros).map(|(w, x)| w * f(x)).sum()
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
