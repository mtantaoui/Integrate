//! Newton-Cotes quadrature rules for numerical integration.
//!
//! Newton-Cotes formulas approximate $\int_a^b f(x)\\,dx$ by replacing $f$ with
//! a polynomial interpolant through equally-spaced nodes and integrating exactly.

use num::{Float, ToPrimitive, Unsigned};

use crate::utils::{
    checkers::check_newton_method_args,
    integrand::{node_count, Integrand, Scalar},
};

/// Converts a node index to `F`.
///
/// Goes through `i64`: on x86-64 a signed conversion is a single instruction,
/// while `usize -> f64` needs a multi-instruction sequence, which is noticeable
/// in the per-node loop. Node counts never approach `i64::MAX`.
#[inline]
fn index<F: Float>(i: usize) -> F {
    F::from(i as i64).unwrap()
}

/// Adds `weight * (f(a) + f(b))` to `sum`, evaluating both endpoints exactly.
///
/// The Newton-Cotes rules sum their interior nodes in one branch-free pass,
/// then the two endpoints.
fn endpoints<F: Float>(f: &mut impl Integrand<F>, a: F, b: F, weight: f64, sum: &mut [f64; 1]) {
    f.accumulate(0..2, |i| if i == 0 { a } else { b }, |_| [weight], sum);
}

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
    rectangle(&mut Scalar(func), lower_limit, upper_limit, n_intervals)
}

/// Composite midpoint rule, shared by [`rectangle_rule`] and its batch variant.
pub(crate) fn rectangle<F, U>(f: &mut impl Integrand<F>, a: F, b: F, n_intervals: U) -> f64
where
    F: Float,
    U: Unsigned + ToPrimitive + Copy,
{
    check_newton_method_args(a, b, n_intervals);

    let n = n_intervals.to_usize().unwrap();
    let h: F = (b - a) / F::from(n).expect("failed to convert n to float");
    let half = F::from(0.5).unwrap();

    let [sum] = f.weighted_sums(n, |i| a + (index::<F>(i) + half) * h, |_| [1.0]);

    sum * h.to_f64().unwrap()
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
    trapezoidal(&mut Scalar(func), lower_limit, upper_limit, n_intervals)
}

/// Composite trapezoidal rule, shared by [`trapezoidal_rule`] and its batch variant.
///
/// Nodes $x_i = a + ih$, $i = 0, \dots, n$ (with $x_n = b$ exactly); weights
/// $\tfrac12, 1, \dots, 1, \tfrac12$.
/// Interior nodes are summed first, then the two endpoints.
pub(crate) fn trapezoidal<F, U>(f: &mut impl Integrand<F>, a: F, b: F, n_intervals: U) -> f64
where
    F: Float,
    U: Unsigned + ToPrimitive + Copy,
{
    check_newton_method_args(a, b, n_intervals);

    let n = n_intervals.to_usize().unwrap();
    let h: F = (b - a) / F::from(n).expect("failed to convert n to float");

    node_count(n, 1); // overflow check for the n + 1 nodes

    let mut sum = [0.0];
    f.accumulate(1..n, |i| a + index::<F>(i) * h, |_| [1.0], &mut sum);
    endpoints(f, a, b, 0.5, &mut sum);

    sum[0] * h.to_f64().unwrap()
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
    simpson(&mut Scalar(f), a, b, n)
}

/// Composite Simpson's 1/3 rule, shared by [`simpson_rule`] and its batch variant.
///
/// Nodes $a + jh/2$, $j = 0, \dots, 2n$ (with the last one $= b$ exactly);
/// weights $1, 4, 2, 4, \dots, 2, 4, 1$, scaled by $h/6$.
/// Interior nodes are summed first, then the two endpoints.
pub(crate) fn simpson<F, U>(f: &mut impl Integrand<F>, a: F, b: F, n: U) -> f64
where
    F: Float,
    U: Unsigned + ToPrimitive + Copy,
{
    check_newton_method_args(a, b, n);

    let n = n.to_usize().unwrap();
    let h: F = (b - a) / F::from(n).expect("failed to convert n to float");
    let h2 = h / F::from(2).unwrap();
    let last = node_count(n, 2) - 1;

    let mut sum = [0.0];
    f.accumulate(
        1..last,
        |j| a + index::<F>(j) * h2,
        |j| if j % 2 == 1 { [4.0] } else { [2.0] },
        &mut sum,
    );
    endpoints(f, a, b, 1.0, &mut sum);

    sum[0] * h.to_f64().unwrap() / 6.0
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
    newton(&mut Scalar(func), lower_limit, upper_limit, n_intervals)
}

/// Composite Newton's 3/8 rule, shared by [`newton_rule`] and its batch variant.
///
/// Nodes $a + jh/3$, $j = 0, \dots, 3n$ (with the last one $= b$ exactly);
/// weights $1, 3, 3, 2, 3, 3, 2, \dots, 3, 3, 1$, scaled by $h/8$.
/// Interior nodes are summed first, then the two endpoints.
pub(crate) fn newton<F, U>(f: &mut impl Integrand<F>, a: F, b: F, n_intervals: U) -> f64
where
    F: Float,
    U: Unsigned + ToPrimitive + Copy,
{
    check_newton_method_args(a, b, n_intervals);

    let n = n_intervals.to_usize().unwrap();
    let h: F = (b - a) / F::from(n).expect("failed to convert n to float");
    let h3 = h / F::from(3).unwrap();
    node_count(n, 3); // overflow check for the 3n + 1 nodes

    // Constant-weight passes (a per-node `j % 3` weight costs ~20%):
    // panel breakpoints a + ih (weight 2), then the two inner nodes of each
    // panel, a + (3i+1)h/3 and a + (3i+2)h/3 (weight 3).
    let mut sum = [0.0];
    f.accumulate(1..n, |i| a + index::<F>(3 * i) * h3, |_| [2.0], &mut sum);
    for offset in [1, 2] {
        f.accumulate(
            0..n,
            |i| a + index::<F>(3 * i + offset) * h3,
            |_| [3.0],
            &mut sum,
        );
    }
    endpoints(f, a, b, 1.0, &mut sum);

    sum[0] * h.to_f64().unwrap() / 8.0
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
    #[test]
    #[should_panic(expected = "Integral limits a and b can't be NaN")]
    fn test_lower_limit_nan() {
        fn f(x: f64) -> f64 {
            x
        }
        let _ = rectangle_rule(f, f64::NAN, 1.0, 1usize);
    }

    #[test]
    #[should_panic(expected = "Integral limits a and b can't be NaN")]
    fn test_upper_limit_nan() {
        fn f(x: f64) -> f64 {
            x
        }
        let _ = rectangle_rule(f, 0.0, f64::NAN, 1usize);
    }
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
    #[test]
    #[should_panic(expected = "Integral limits a and b can't be NaN")]
    fn test_lower_limit_nan() {
        fn f(x: f64) -> f64 {
            x
        }
        let _ = trapezoidal_rule(f, f64::NAN, 1.0, 1usize);
    }

    #[test]
    #[should_panic(expected = "Integral limits a and b can't be NaN")]
    fn test_upper_limit_nan() {
        fn f(x: f64) -> f64 {
            x
        }
        let _ = trapezoidal_rule(f, 0.0, f64::NAN, 1usize);
    }
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

    #[test]
    #[should_panic(expected = "Integral limits a and b can't be NaN")]
    fn test_lower_limit_nan() {
        fn f(x: f64) -> f64 {
            x
        }
        let _ = simpson_rule(f, f64::NAN, 1.0, 1usize);
    }

    #[test]
    #[should_panic(expected = "Integral limits a and b can't be NaN")]
    fn test_upper_limit_nan() {
        fn f(x: f64) -> f64 {
            x
        }
        let _ = simpson_rule(f, 0.0, f64::NAN, 1usize);
    }
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

    #[test]
    #[should_panic(expected = "Integral limits a and b can't be NaN")]
    fn test_lower_limit_nan() {
        fn f(x: f64) -> f64 {
            x
        }
        let _ = newton_rule(f, f64::NAN, 1.0, 1usize);
    }

    #[test]
    #[should_panic(expected = "Integral limits a and b can't be NaN")]
    fn test_upper_limit_nan() {
        fn f(x: f64) -> f64 {
            x
        }
        let _ = newton_rule(f, 0.0, f64::NAN, 1usize);
    }
}
