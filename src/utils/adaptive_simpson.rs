use std::{fmt, ops::MulAssign};

use num::Float;

/// Represents a subinterval $[a, b]$ during adaptive subdivision, storing the function
/// values at five evenly-spaced points and a linked pointer to the next pending subinterval.
#[derive(Clone, Debug)]
pub struct SubInterval<F: Float> {
    /// Upper limit of the subinterval $b$.
    pub upper_limit: F,
    /// Lower limit of the subinterval $a$.
    pub lower_limit: F,
    /// Function evaluations at the five evenly-spaced points:
    /// $a$, $a+h/4$, $(a+b)/2$, $b-h/4$, $b$,
    /// where $h = b - a$.
    pub function: [F; 5],
    /// Next pending subinterval in the stack; `None` if this is the last one.
    pub interval: Option<Box<SubInterval<F>>>,
}

/// Error returned when adaptive Simpson fails to meet tolerance before the minimum step size
/// `min_h` is reached.
#[derive(Debug, Clone)]
pub struct AdaptiveSimpsonError;

impl fmt::Display for AdaptiveSimpsonError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let msg = "No subinterval of length > min_h was found for which the estimated error was less that the pro-rated tolerance";
        write!(f, "{}", msg)
    }
}

/// Given a subinterval, evaluates the two missing function points and returns `(s1, s2)`
/// where `s1` is the single-interval Simpson estimate and `s2` is the composite
/// (two sub-subintervals) Simpson estimate.  The ratio $|s_1 - s_2|$ is used as a
/// local error indicator.
pub fn simpson_rule_update<Func, F: Float + MulAssign + fmt::Debug>(
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
