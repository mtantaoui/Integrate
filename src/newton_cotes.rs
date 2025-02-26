use std::ops::Div;

use num::{Float, ToPrimitive, Unsigned};
use rayon::prelude::*;

use crate::utils::checkers::check_newton_method_args;

/// This function integrates $f(x)$ from $a$ to $a+nh$ using the rectangle
/// rule by summing from the left end of the interval to the right end.
///
/// * `func` - Integrand function of a single variable.
/// * `lower_limit` - lower limit of the integration interval.
/// * `upper_limit` - upper limit of the integration interval.
/// * `n_intervals` - number of subintervals.
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
pub fn rectangle_rule<Func, F1: Float + Sync, F2: Float + Sync, U: Unsigned + ToPrimitive + Copy>(
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

    let integral: f64 = (0..(n_intervals.to_usize().unwrap()))
        .into_par_iter()
        .map(|i| {
            // subinterval index (as real)
            let i = F1::from(i).expect("failed to convert subinterval index i");

            // subinterval midpoint
            let x = lower_limit
                + i * h
                + (h / F1::from(2).expect("failed to compute subinterval midpoint"));

            // converting f(x) to primitive type f64
            func(x).to_f64().expect("failed to convert f(x) to f64")
            // 0.0
        })
        .sum();
    integral * h.to_f64().unwrap()
}

/// This function integrates $f(x)$ from $a$ to $a+nh$ using the Simpson's
/// rule by summing from the left end of the interval to the right end.
///  
/// * `func` - Integrand function of a single variable.
/// * `lower_limit` - lower limit of the integration interval.
/// * `upper_limit` - upper limit of the integration interval.
/// * `n_intervals` - number of subintervals.
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
pub fn trapezoidal_rule<
    Func,
    F1: Float + Sync,
    F2: Float + Send,
    U: Unsigned + ToPrimitive + Copy,
>(
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

    // first term of the sum
    let i_0 = func(lower_limit).to_f64().unwrap();

    let integral: f64 = (1..(n_intervals.to_usize().unwrap()))
        .into_par_iter()
        .map(|i| {
            // subinterval index (as real)
            let i = F1::from(i).expect("failed to convert subinterval index i");
            func(lower_limit + i * h).to_f64().unwrap()
        })
        .sum();

    let n: F1 = F1::from(n_intervals).expect("failed to convert number of steps n");
    // last term of the sum
    let i_n = func(lower_limit + h * n).to_f64().unwrap();

    (0.5 * i_0 + integral + 0.5 * i_n) * h.to_f64().expect("failed to convert subintervql length")
}

/// This function integrates $f(x)$ from $a$ to $a+nh$ using the Simpson's
/// rule by summing from the left end of the interval to the right end.
///
/// * `func` - Integrand function of a single variable.
/// * `lower_limit` - lower limit of the integration interval.
/// * `upper_limit` - upper limit of the integration interval.
/// * `n_intervals` - number of subintervals.
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
pub fn simpson_rule<Func, F1: Float + Sync, F2: Float, U: Unsigned + ToPrimitive + Copy>(
    f: Func,
    a: F1,
    b: F1,
    n: U,
) -> f64
where
    Func: Fn(F1) -> F2 + Sync,
{
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
mod rectangle_rule_tests {
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

        let square = |x: f32| x * x;

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

#[cfg(test)]
mod trapezoidal_rule_tests {
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
