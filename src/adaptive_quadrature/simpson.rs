//! Adaptive Quadrature
//!
//! If an integrand is poorly behaved in a small interval about a point,
//! then an attempt to integrate the function over an interval which contains
//! the poorly behaved interval either requires that small subintervals
//! are chosen for composite quadratures or the interval is decomposed into three intervals,
//! two on which the function is well-behaved and relatively large subintervals
//! can be chosen for the composite quadrature technique and one in which smaller subintervals need to be chosen.
//!
//! Adaptive techniques are attempts to automatically detect and control the length of subintervals.
//!
//! The technique for which the link to the listing is given below uses Simpson's rule
//! for integrating a function $f(x)$ on a closed and bounded interval $\[a,b\]$.

use num::Float;
use std::fmt;

use std::ops::{AddAssign, MulAssign};

#[derive(Clone, Debug)]
struct SubInterval<F: Float> {
    upper_limit: F,
    lower_limit: F,
    function: [F; 5],
    interval: Option<Box<SubInterval<F>>>,
}

type Result<T> = std::result::Result<T, AdaptiveSimpsonError>;

#[derive(Debug, Clone)]
pub struct AdaptiveSimpsonError;

impl fmt::Display for AdaptiveSimpsonError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let msg = "No subinterval of length > min_h was found for which the estimated error was less that the pro-rated tolerance";
        write!(f, "{}", msg)
    }
}
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
/// use integrate::adaptive_quadrature::simpson::adaptive_simpson_method;
///
///
/// fn f(x: f64) -> f64 {
///     x.exp()
/// }
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

fn simpson_rule_update<Func, F: Float + MulAssign + fmt::Debug>(
    func: Func,
    pinterval: &mut SubInterval<F>,
) -> (F, F)
where
    Func: Fn(F) -> F + Sync,
{
    let two = F::one() + F::one();
    let four = two + two;
    let six = four + two;

    let h = pinterval.upper_limit - pinterval.lower_limit;
    let h4 = h / four;

    pinterval.function[1] = func(pinterval.lower_limit + h4);
    pinterval.function[3] = func(pinterval.upper_limit - h4);

    let mut s1 = pinterval.function[0] + four * pinterval.function[2] + pinterval.function[4];
    s1 *= h / six;

    let mut s2 = pinterval.function[0]
        + four * pinterval.function[1]
        + two * pinterval.function[2]
        + four * pinterval.function[3]
        + pinterval.function[4];
    s2 *= h / (six * two);

    (s1, s2)
}

// tests in tests/test_adaptive_quadrature.rs
