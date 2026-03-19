//! Newton-Cotes quadrature rules for numerical integration.
//!
//! Newton-Cotes formulas approximate $\int_a^b f(x)\\,dx$ by replacing $f$ with
//! a polynomial interpolant through equally-spaced nodes and integrating exactly.

use num::{Float, ToPrimitive, Unsigned};

use crate::utils::checkers::check_newton_method_args;

/// Approximates $\int_a^b f(x)\\,dx$ using the composite midpoint rectangle rule.
///
/// Each subinterval $[a + ih,\\; a + (i+1)h]$ contributes $h \cdot f$ evaluated at
/// its midpoint $a + (i + \tfrac{1}{2})h$.  The composite formula is
///
/// $$R_h(f) = h \sum_{i=0}^{n-1} f\left(a + \left(i + \tfrac{1}{2}\right)h\right), \quad h = \frac{b-a}{n}.$$
///
/// **Truncation error** (requires $f \in C^2\[a,b\]$): for some $c \in \[a,b\]$,
///
/// $$R_h(f) - \int_a^b f(x)\\,dx = -\frac{h^2}{24}(b-a)f^{\prime\prime}(c).$$
///
/// # Parameters
/// * `func` - Integrand function of a single variable.
/// * `lower_limit` - Lower limit $a$ of the integration interval.
/// * `upper_limit` - Upper limit $b$ of the integration interval.
/// * `n_intervals` - Number of subintervals $n$.
///
/// # Panics
/// Delegates argument validation to `check_newton_method_args`,
/// which panics if `n_intervals` is zero, either limit is infinite or NaN, or `lower_limit > upper_limit`.
///
/// # Examples
/// ```
/// use integrate::newton_cotes::rectangle_rule;
///
///
/// let square = |x: f64| x * x;
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
pub fn rectangle_rule<Func, F1: Float, F2: Float, U: Unsigned + ToPrimitive + Copy>(
    func: Func,
    lower_limit: F1,
    upper_limit: F1,
    n_intervals: U,
) -> f64
where
    Func: Fn(F1) -> F2,
{
    check_newton_method_args(lower_limit, upper_limit, n_intervals);

    let n = n_intervals.to_usize().unwrap();
    let h: F1 = (upper_limit - lower_limit) / F1::from(n).expect("failed to convert n to float");
    let half = F1::from(0.5).unwrap();

    let integral: f64 = (0..n)
        .map(|i| {
            let x = lower_limit + (F1::from(i).unwrap() + half) * h;
            func(x).to_f64().expect("failed to convert f(x) to f64")
        })
        .sum();

    integral * h.to_f64().unwrap()
}

/// Approximates $\int_a^b f(x)\\,dx$ using the composite trapezoidal rule.
///
/// Each subinterval is approximated by a trapezoid.  The composite formula is
///
/// $$T_h(f) = h\left[\tfrac{1}{2}f(a) + f(a+h) + f(a+2h) + \cdots + f(b-h) + \tfrac{1}{2}f(b)\right], \quad h = \frac{b-a}{n}.$$
///
/// **Truncation error** (requires $f \in C^2\[a,b\]$): for some $c \in \[a,b\]$,
///
/// $$T_h(f) - \int_a^b f(x)\\,dx = -\frac{h^2}{12}(b-a)f^{\prime\prime}(c).$$
///
/// # Parameters
/// * `func` - Integrand function of a single variable.
/// * `lower_limit` - Lower limit $a$ of the integration interval.
/// * `upper_limit` - Upper limit $b$ of the integration interval.
/// * `n_intervals` - Number of subintervals $n$.
///
/// # Panics
/// Delegates argument validation to `check_newton_method_args`,
/// which panics if `n_intervals` is zero, either limit is infinite or NaN, or `lower_limit > upper_limit`.
///
/// # Examples
/// ```
/// use integrate::newton_cotes::trapezoidal_rule;
///
///
/// let square = |x: f64| x * x;
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
pub fn trapezoidal_rule<Func, F1: Float, F2: Float, U: Unsigned + ToPrimitive + Copy>(
    func: Func,
    lower_limit: F1,
    upper_limit: F1,
    n_intervals: U,
) -> f64
where
    Func: Fn(F1) -> F2,
{
    check_newton_method_args(lower_limit, upper_limit, n_intervals);

    let n = n_intervals.to_usize().unwrap();
    let h: F1 = (upper_limit - lower_limit) / F1::from(n).expect("failed to convert n to float");

    let interior: f64 = (1..n)
        .map(|i| {
            func(lower_limit + F1::from(i).unwrap() * h)
                .to_f64()
                .unwrap()
        })
        .sum();

    let endpoints =
        0.5 * (func(lower_limit).to_f64().unwrap() + func(upper_limit).to_f64().unwrap());

    (endpoints + interior) * h.to_f64().unwrap()
}

/// Approximates $\int_a^b f(x)\\,dx$ using the composite Simpson's rule.
///
/// Each pair of subintervals is approximated by a quadratic (parabolic) arc.
/// The composite formula is
///
/// $$S_h(f) = \frac{h}{6}\left[f(a) + 4f\left(a+\tfrac{h}{2}\right) + 2f(a+h) + 4f\left(a+\tfrac{3h}{2}\right) + \cdots + f(b)\right], \quad h = \frac{b-a}{n}.$$
///
/// **Truncation error** (requires $f \in C^4\[a,b\]$): for some $c \in \[a,b\]$,
///
/// $$S_h(f) - \int_a^b f(x)\\,dx = -\frac{h^4}{180}(b-a)f^{(4)}(c).$$
///
/// Simpson's rule is exact for polynomials of degree $\leq 3$.
///
/// # Parameters
/// * `func` - Integrand function of a single variable.
/// * `lower_limit` - Lower limit $a$ of the integration interval.
/// * `upper_limit` - Upper limit $b$ of the integration interval.
/// * `n_intervals` - Number of subintervals $n$.
///
/// # Panics
/// Delegates argument validation to `check_newton_method_args`,
/// which panics if `n_intervals` is zero, either limit is infinite or NaN, or `lower_limit > upper_limit`.
///
/// # Examples
/// ```
/// use integrate::newton_cotes::simpson_rule;
///
///
/// let square = |x: f64| x * x;
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
pub fn simpson_rule<Func, F1: Float, F2: Float, U: Unsigned + ToPrimitive + Copy>(
    f: Func,
    a: F1,
    b: F1,
    n: U,
) -> f64
where
    Func: Fn(F1) -> F2,
{
    check_newton_method_args(a, b, n);

    let n_usize = n.to_usize().unwrap();
    let h: F1 = (b - a) / F1::from(n_usize).unwrap();
    let h2 = h / F1::from(2).unwrap();

    // f(a) + 4f(a + h/2)
    let first = f(a).to_f64().unwrap() + 4.0 * f(a + h2).to_f64().unwrap();

    // Σ 2f(a+ih) + 4f(a+ih+h/2)  for i = 1..n-1, then 2f(a+(n-1)h) + 4f(a+(n-1)h+h/2)
    let middle: f64 = (1..n_usize)
        .map(|i| {
            let xi = a + F1::from(i).unwrap() * h;
            2.0 * f(xi).to_f64().unwrap() + 4.0 * f(xi + h2).to_f64().unwrap()
        })
        .sum();

    let last = f(b).to_f64().unwrap();

    (first + middle + last) * h.to_f64().unwrap() / 6.0
}

/// Approximates $\int_a^b f(x)\\,dx$ using the composite Newton's 3/8 rule.
///
/// Each subinterval uses four equally-spaced points (a cubic interpolant).
/// The composite formula is
///
/// $$N_h(f) = \frac{3h}{8}\left[f(a) + 3f\left(a+\tfrac{h}{3}\right) + 3f\left(a+\tfrac{2h}{3}\right) + 2f(a+h) + \cdots + f(b)\right], \quad h = \frac{b-a}{n}.$$
///
/// **Truncation error** (requires $f \in C^4\[a,b\]$): for some $c \in \[a,b\]$,
///
/// $$N_h(f) - \int_a^b f(x)\\,dx = -\frac{3h^4}{80}(b-a)f^{(4)}(c).$$
///
/// Newton's 3/8 rule is exact for polynomials of degree $\leq 3$.
///
/// # Parameters
/// * `func` - Integrand function of a single variable.
/// * `lower_limit` - Lower limit $a$ of the integration interval.
/// * `upper_limit` - Upper limit $b$ of the integration interval.
/// * `n_intervals` - Number of subintervals $n$.
///
/// # Panics
/// Delegates argument validation to `check_newton_method_args`,
/// which panics if `n_intervals` is zero, either limit is infinite or NaN, or `lower_limit > upper_limit`.
///
/// # Examples
/// ```
/// use integrate::newton_cotes::newton_rule;
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
pub fn newton_rule<Func, F1: Float, F2: Float, U: Unsigned + ToPrimitive + Copy>(
    func: Func,
    lower_limit: F1,
    upper_limit: F1,
    n_intervals: U,
) -> f64
where
    Func: Fn(F1) -> F2,
{
    check_newton_method_args(lower_limit, upper_limit, n_intervals);

    let n = n_intervals.to_usize().unwrap();
    let h: F1 = (upper_limit - lower_limit) / F1::from(n).unwrap();
    let h3 = h / F1::from(3).unwrap();

    // f(a) + 3f(a+h/3) + 3f(a+2h/3)
    let first = func(lower_limit).to_f64().unwrap()
        + 3.0
            * (func(lower_limit + h3).to_f64().unwrap()
                + func(lower_limit + F1::from(2).unwrap() * h3)
                    .to_f64()
                    .unwrap());

    // For i = 1..n-1: 2f(a+ih) + 3f(a+ih+h/3) + 3f(a+ih+2h/3)
    let middle: f64 = (1..n)
        .map(|i| {
            let xi = lower_limit + F1::from(i).unwrap() * h;
            2.0 * func(xi).to_f64().unwrap()
                + 3.0
                    * (func(xi + h3).to_f64().unwrap()
                        + func(xi + F1::from(2).unwrap() * h3).to_f64().unwrap())
        })
        .sum();

    let last = func(upper_limit).to_f64().unwrap();

    (first + middle + last) * h.to_f64().unwrap() / 8.0
}

#[cfg(test)]
mod rectangle_rule_tests {
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

        let analytic_result: f64 = 1.0 / 3.0;

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

        let analytic_result: f64 = 1.0 / 3.0;

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

        let analytic_result: f64 = 1.0 / 3.0;

        assert!((integral - analytic_result).abs() < EPSILON);
    }

    #[test]
    fn test_f32_to_f32() {
        // f32 to f32

        let square = |x: f32| x * x;

        let a = 0.0;
        let b = 1.0;

        let integral = rectangle_rule(square, a, b, NUM_STEPS);

        let analytic_result: f64 = 1.0 / 3.0;

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

#[cfg(test)]
mod trapezoidal_rule_tests {
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

        let analytic_result: f64 = 1.0 / 3.0;

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

        let analytic_result: f64 = 1.0 / 3.0;

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

        let analytic_result: f64 = 1.0 / 3.0;

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

        let analytic_result: f64 = 1.0 / 3.0;

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

#[cfg(test)]
mod simpson_rule_tests {

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

        let analytic_result: f64 = 1.0 / 3.0;

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

        let analytic_result: f64 = 1.0 / 3.0;

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

        let analytic_result: f64 = 1.0 / 3.0;

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

        let analytic_result: f64 = 1.0 / 3.0;

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

#[cfg(test)]
mod newton_rule_tests {

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

        let analytic_result: f64 = 1.0 / 3.0;

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

        let analytic_result: f64 = 1.0 / 3.0;

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

        let analytic_result: f64 = 1.0 / 3.0;

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

        let analytic_result: f64 = 1.0 / 3.0;

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
