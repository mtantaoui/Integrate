//! Adaptive Simpson quadrature for numerical integration.
//!
//! Unlike fixed-step methods, adaptive quadrature automatically refines the
//! step size in regions where the integrand varies rapidly, achieving a
//! user-specified tolerance with fewer function evaluations.

use std::{
    fmt,
    ops::{AddAssign, MulAssign},
};

use num::Float;

use crate::utils::adaptive_simpson::{simpson_rule_update, AdaptiveSimpsonError, SubInterval};

type Result<T> = std::result::Result<T, AdaptiveSimpsonError>;

/// Numerically integrates $f$ over $[a, b]$ using the adaptive Simpson method.
///
/// Starting at the left endpoint $a$, the algorithm finds the minimum power of 2, $m$,
/// such that the difference between the single-interval Simpson estimate and the
/// composite (two sub-subintervals) Simpson estimate on $\left[a,\\, a+\frac{b-a}{2^m}\right]$
/// is less than the pro-rated tolerance:
///
/// $$2 \cdot \texttt{tolerance} \cdot \frac{\text{length of subinterval}}{b - a}$$
///
/// The process then repeats for the remaining interval $\left[a + \frac{b-a}{2^m},\\, b\right]$,
/// advancing until the right endpoint $b$ is reached.  The integral is the sum of the
/// accepted subinterval estimates.
///
/// # Parameters
///
///  * `func` — Integrand $f$; a single-variable function.
///  * `lower_limit` — Lower limit of integration $a$.
///  * `upper_limit` — Upper limit of integration $b$; must satisfy `upper_limit > lower_limit`.
///  * `min_h` — Minimum allowed subinterval length; must be a positive finite float.
///  * `tolerance` — Desired absolute error bound; must be a positive finite float.
///
/// # Errors
///
/// Returns `Err(AdaptiveSimpsonError)` if a subinterval of length `min_h` is reached before the
/// tolerance is satisfied. This typically means the integrand has a singularity or varies too
/// rapidly — try increasing `tolerance` or decreasing `min_h`.
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
    assert!(
        min_h > F::zero() && min_h.is_finite(),
        "min_h must be positive and finite"
    );
    assert!(
        tolerance > F::zero() && tolerance.is_finite(),
        "tolerance must be positive and finite"
    );
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
