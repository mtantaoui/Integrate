use std::{
    fmt,
    ops::{AddAssign, MulAssign},
};

use num::Float;

use crate::utils::adaptive_simpson::{simpson_rule_update, AdaptiveSimpsonError, SubInterval};

type Result<T> = std::result::Result<T, AdaptiveSimpsonError>;

/// Simpson-Simpson adaptive method
///
/// Integrate, using the Simpson-Simpson adaptive method, the user supplied function $f$ from $a$ to $b$.
///
/// * `func` - Integrand function of a single variable.
/// * `lower_limit` is the lower limit of integration.
/// * `upper_limit`  is the upper limit of integration where `upper_limit` > `lower_limit`.
/// * `tolerance` is the tolerance.
/// * `min_h` is the minimum subinterval length to be used.
///
/// Starting at the left-end point, `a`, find the min power of 2, $m$,
/// so that the difference between using Simpson's rule and the composite Simpson's
/// rule on the interval $\[a, a+\frac{b-a}{2^m}\]$ is less than
/// $$ 2 * \verb|tolerance| * \frac{\text{length of the subinterval}}{b-a}$$
///
/// Then repeat the process for integrating over the interval $\[a + \frac{b-a}{2^m}, b\]$ until the right
/// end point, b, is finally reached.
///
/// The integral is then the sum of the integrals of each subinterval.  If at any time,
/// the length of the subinterval for which the estimates based on Simpson's rule and
/// the composite Simpson's rule is less than `min_h`, the process is terminated with an
/// `AdaptiveSimpsonError` error.
///
/// # Examples
/// ```
/// use integrate::adaptive_quadrature::adaptive_simpson_method;
///
///
/// let f = |x: f64| x.exp();
///
/// let a = 0.0;
/// let b = 1.0;
///
/// let tolerance = 10.0e-6;
/// let min_h = 10.0e-3;
///
///
/// let result = adaptive_simpson_method(f, a, b, min_h, tolerance);
///
///
/// match result{
///     Ok(res)=>{println!("{}", res)}
///     Err(err)=>{println!("{}", err)}
/// };
///
/// ```
pub fn adaptive_simpson_method<Func, F: Float + MulAssign + AddAssign + fmt::Debug>(
    func: Func,
    lower_limit: F,
    upper_limit: F,
    min_h: F,
    tolerance: F,
) -> Result<F>
where
    Func: Fn(F) -> F + Sync + Copy,
{
    let two = F::one() + F::one();

    let mut integral: F = F::zero();
    let epsilon_density = two * tolerance / (upper_limit - lower_limit);

    // Create the initial level, with lower_limit = a, upper_limit = b,
    // and f(x) evaluated at a, b, and (a + b) / 2.

    let interval: SubInterval<F> = SubInterval {
        upper_limit,
        lower_limit,
        function: [
            func(lower_limit),
            F::nan(),
            func((lower_limit + upper_limit) / two),
            F::nan(),
            func(upper_limit),
        ],
        interval: None,
    };

    let mut pinterval = Box::new(interval);

    // Calculate the tolerance for the current interval.
    // calculate the single subinterval Simpson rule,
    // and the two subintervals composite Simpson rule.

    let mut epsilon = epsilon_density * (upper_limit - lower_limit);
    let (mut s1, mut s2) = simpson_rule_update(func, &mut pinterval);

    let mut qinterval: SubInterval<F>;

    while pinterval.upper_limit - pinterval.lower_limit > min_h {
        if (s1 - s2).abs() < epsilon {
            // If the two estimates are close, then increment the
            // integral and if we are not at the right end, set the
            // left end of the new interval to the right end of the
            // old interval and the right end of the new interval
            // remains the same (as the previous right end for this
            // interval.

            integral += s2;

            if pinterval.interval.is_none() {
                return Ok(integral);
            }

            // Move to the next interval
            qinterval = *pinterval.interval.take().unwrap();
            qinterval.lower_limit = pinterval.upper_limit;
            qinterval.function[0] = qinterval.function[2];
            qinterval.function[2] = qinterval.function[3];

            pinterval = Box::new(qinterval);
        } else {
            // If the two estimates are not close, then create a new
            // interval with same left end point and right end point
            // at the midpoint of the current interval.

            let limit1 = pinterval.lower_limit;
            let limit2 = (pinterval.upper_limit + pinterval.lower_limit) / two;

            let upper_limit = if limit1 > limit2 { limit1 } else { limit2 };
            let lower_limit = if limit1 > limit2 { limit2 } else { limit1 };

            qinterval = SubInterval {
                lower_limit,
                upper_limit,
                function: [F::nan(); 5],
                interval: None,
            };

            qinterval.function[0] = pinterval.function[0];
            qinterval.function[2] = pinterval.function[1];
            qinterval.function[4] = pinterval.function[2];

            qinterval.interval = Some(pinterval);

            pinterval = Box::new(qinterval);
        }

        // Update Simpson's rule for the new interval
        (s1, s2) = simpson_rule_update(func, &mut pinterval);
        epsilon = epsilon_density * (pinterval.upper_limit - pinterval.lower_limit);
    }
    Err(AdaptiveSimpsonError)
}

// tests in tests/test_adaptive_quadrature.rs
