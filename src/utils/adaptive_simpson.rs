use std::{fmt, ops::MulAssign};

use num::Float;

#[derive(Clone, Debug)]
pub struct SubInterval<F: Float> {
    pub upper_limit: F,
    pub lower_limit: F,
    pub function: [F; 5],
    pub interval: Option<Box<SubInterval<F>>>,
}

#[derive(Debug, Clone)]
pub struct AdaptiveSimpsonError;

impl fmt::Display for AdaptiveSimpsonError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let msg = "No subinterval of length > min_h was found for which the estimated error was less that the pro-rated tolerance";
        write!(f, "{}", msg)
    }
}

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
