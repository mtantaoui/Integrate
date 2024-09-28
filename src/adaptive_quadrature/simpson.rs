use num::{Float, ToPrimitive, Unsigned};
use std::fmt::Debug;
use std::ops::{AddAssign, MulAssign};

#[derive(Clone, Debug)]
struct SubInterval<F: Float> {
    upper_limit: F,
    lower_limit: F,
    function: [F; 5],
    interval: Option<Box<SubInterval<F>>>,
}

pub fn adaptive_simpson_method<
    F: Float + MulAssign + AddAssign + Debug,
    U: Unsigned + ToPrimitive + Copy,
>(
    f: fn(F) -> F,
    a: F,
    b: F,
    n: U,
) -> F {
    let two = F::one() + F::one();
    let tolerance = F::from(f32::EPSILON).unwrap();

    let mut integral: F = F::zero();
    let epsilon_density = two * tolerance / (b - a);

    // Create the initial level, with lower_limit = a, upper_limit = b,
    // and f(x) evaluated at a, b, and (a + b) / 2.

    let interval: SubInterval<F> = SubInterval {
        upper_limit: b,
        lower_limit: a,
        function: [f(a), F::nan(), f((a + b) / two), F::nan(), f(b)],
        interval: None,
    };

    let mut pinterval = Box::new(interval);

    // Calculate the tolerance for the current interval.
    // calculate the single subinterval Simpson rule,
    // and the two subintervals composite Simpson rule.

    let mut epsilon = epsilon_density * (b - a);
    let (mut s1, mut s2) = simpson_rule_update(f, &mut pinterval);

    let n = F::from(n).unwrap();
    let min_h = (b - a) / n;

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
                return integral;
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
            qinterval = SubInterval {
                lower_limit: pinterval.lower_limit,
                upper_limit: (pinterval.upper_limit + pinterval.lower_limit) / two,
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
        (s1, s2) = simpson_rule_update(f, &mut pinterval);
        epsilon = epsilon_density * (pinterval.upper_limit - pinterval.lower_limit);
    }

    integral
}

fn simpson_rule_update<F: Float + MulAssign + Debug>(
    f: fn(F) -> F,
    pinterval: &mut SubInterval<F>,
) -> (F, F) {
    let two = F::one() + F::one();
    let four = two + two;
    let six = four + two;

    let h = pinterval.upper_limit - pinterval.lower_limit;
    let h4 = h / four;

    pinterval.function[1] = f(pinterval.lower_limit + h4);
    pinterval.function[3] = f(pinterval.upper_limit - h4);

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
