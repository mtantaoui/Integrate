use std::{f64::consts::PI, iter::Sum};

use num::Float;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

const EPSILON: f64 = 10e-4;
const NUM_STEPS: usize = 100_000;

#[derive(Clone, Debug)]
pub struct Problem<F: Float> {
    pub id: usize,
    pub function: fn(F) -> F, // represents integrand function
    pub limits: (F, F),       // integration limits
    pub exact: F,             // exact value of the integral
    pub n: usize,             // usually the number of steps
}

impl<F: Float> Problem<F> {
    pub fn check_result(&self, result: F) -> bool {
        let epsilon = F::from(EPSILON).unwrap();
        (self.exact - result).abs() < epsilon
    }
}

pub fn problem01<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        x.exp()
    }

    Problem {
        id: 1,
        function: f,
        limits: (F::zero(), F::one()),
        exact: F::one().exp() - F::one(),
        n: NUM_STEPS,
    }
}

pub fn problem02<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        let atom = F::from(0.3).unwrap();

        if x < atom {
            F::zero()
        } else {
            F::one()
        }
    }

    let exact = F::from(0.7).unwrap();

    Problem {
        id: 2,
        function: f,
        limits: (F::zero(), F::one()),
        exact,
        n: NUM_STEPS,
    }
}

pub fn problem03<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        x.sqrt()
    }

    let exact = F::from(2.0 / 3.0).unwrap();

    Problem {
        id: 3,
        function: f,
        limits: (F::zero(), F::one()),
        exact,
        n: NUM_STEPS,
    }
}

pub fn problem04<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        let constant = F::from(0.92).unwrap();
        constant * x.cosh() - x.cos()
    }

    let exact = F::from(1.84 * 1.0.sinh() - 2.0 * 1.0.sin()).unwrap();

    Problem {
        id: 4,
        function: f,
        limits: (-F::one(), F::one()),
        exact,
        n: NUM_STEPS,
    }
}

pub fn problem05<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        let constant = F::from(0.9).unwrap();
        F::one() / (x.powi(4) + x.powi(2) + constant)
    }

    let exact = F::from(1.582_232_963_729_673).unwrap();

    Problem {
        id: 5,
        function: f,
        limits: (-F::one(), F::one()),
        exact,
        n: NUM_STEPS,
    }
}

pub fn problem06<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        let constant = F::from(0.5).unwrap();
        (x + constant).abs().sqrt()
    }

    let exact = (2.0.sqrt() + 3.0 * 6.0.sqrt()) / 6.0;
    let exact = F::from(exact).unwrap();

    Problem {
        id: 6,
        function: f,
        limits: (-F::one(), F::one()),
        exact,
        n: NUM_STEPS,
    }
}

pub fn problem07<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        F::one() / x.sqrt()
    }

    let constant = F::from(10_000).unwrap();
    let two = F::one() + F::one();

    let exact = F::one() + F::one() - two * (F::one() / constant).sqrt();

    Problem {
        id: 7,
        function: f,
        limits: (F::one() / constant, F::one()),
        exact,
        n: NUM_STEPS,
    }
}

pub fn problem08<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        F::one() / (F::one() + x.powi(4))
    }

    let exact = F::from(0.866_972_987_339_911).unwrap();

    Problem {
        id: 8,
        function: f,
        limits: (F::zero(), F::one()),
        exact,
        n: NUM_STEPS,
    }
}

pub fn problem09<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        let two = F::one() + F::one();
        let ten = F::from(10).unwrap();
        let pi = F::from(PI).unwrap();

        two / (two + (ten * pi * x).sin())
    }

    let exact = F::from(2.0 / 3.0.sqrt()).unwrap();

    Problem {
        id: 9,
        function: f,
        limits: (F::zero(), F::one()),
        exact,
        n: NUM_STEPS,
    }
}

pub fn problem10<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        F::one() / (F::one() + x)
    }

    let exact = F::from(2.0.ln()).unwrap();

    Problem {
        id: 10,
        function: f,
        limits: (F::zero(), F::one()),
        exact,
        n: NUM_STEPS,
    }
}

pub fn problem11<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        F::one() / (F::one() + x.exp())
    }

    let exact = F::from(((2.0 * 1.0.exp()) / (1.0 + 1.0.exp())).ln()).unwrap();

    Problem {
        id: 11,
        function: f,
        limits: (F::zero(), F::one()),
        exact,
        n: NUM_STEPS,
    }
}

pub fn problem12<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        x / (x.exp() - F::one())
    }

    let exact = F::from(0.777_504_634_112_248_2).unwrap();

    Problem {
        id: 12,
        function: f,
        limits: (F::epsilon(), F::one()),
        exact,
        n: NUM_STEPS,
    }
}

#[allow(dead_code)]
pub fn problem13<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        x.sin() / x
    }

    let exact = F::from(1.658_347_594_218_874_1).unwrap();

    let ten = F::from(10).unwrap();

    Problem {
        id: 13,
        function: f,
        limits: (F::epsilon(), ten),
        exact,
        n: NUM_STEPS,
    }
}

pub fn problem14<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        let constant = F::from(50).unwrap();
        let pi = F::from(PI).unwrap();
        constant.sqrt() * (-constant * pi * x * x).exp()
    }

    let exact = F::from(0.500000211166).unwrap();

    let ten = F::from(10).unwrap();

    Problem {
        id: 14,
        function: f,
        limits: (F::zero(), ten),
        exact,
        n: NUM_STEPS,
    }
}

pub fn problem15<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        let constant = F::from(25).unwrap();
        constant * (-constant * x).exp()
    }

    let exact = F::from(1.0 - (-250.0).exp()).unwrap();

    let ten = F::from(10).unwrap();

    Problem {
        id: 15,
        function: f,
        limits: (F::zero(), ten),
        exact,
        n: NUM_STEPS,
    }
}

pub fn problem16<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        let constant1 = F::from(50).unwrap();
        let constant2 = F::from(2500).unwrap();

        let pi = F::from(PI).unwrap();

        constant1 / (pi * (constant2 * x * x + F::one()))
    }

    let exact = F::from(0.499_363_381_076_456_7).unwrap();

    let ten = F::from(10).unwrap();

    Problem {
        id: 16,
        function: f,
        limits: (F::zero(), ten),
        exact,
        n: NUM_STEPS,
    }
}

pub fn problem17<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        let constant = F::from(50).unwrap();

        let pi = F::from(PI).unwrap();

        (constant * pi * x).sin().powi(2)
    }

    let exact = F::from(0.5).unwrap();

    Problem {
        id: 17,
        function: f,
        limits: (F::zero(), F::one()),
        exact,
        n: NUM_STEPS,
    }
}

pub fn problem18<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        x / (x.exp() + F::one())
    }

    let exact = F::from(0.170_557_349_502_438_2).unwrap();

    Problem {
        id: 18,
        function: f,
        limits: (F::zero(), F::one()),
        exact,
        n: NUM_STEPS,
    }
}

pub fn problem19<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        x.ln()
    }

    fn anti_derivative_f<F: Float>(x: F) -> F {
        x * x.ln() - x
    }

    let constant = F::one() / F::from(10_000).unwrap();
    let exact = anti_derivative_f(F::one()) - anti_derivative_f(constant);

    Problem {
        id: 19,
        function: f,
        limits: (constant, F::one()),
        exact,
        n: NUM_STEPS,
    }
}

pub fn problem20<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        let constant = F::from(1.005).unwrap();
        F::one() / (x.powi(2) + constant)
    }

    let exact = F::from(1.564_396_444_069_049_9).unwrap();

    Problem {
        id: 20,
        function: f,
        limits: (-F::one(), F::one()),
        exact,
        n: NUM_STEPS,
    }
}

pub fn problem21<F: Float>() -> Problem<F> {
    // evaluates the hyperbolic secant, while avoiding COSH overflow.
    fn sech<F: Float>(x: F) -> F {
        let log_huge = F::from(80.0).unwrap();

        if log_huge < x.abs() {
            F::zero()
        } else {
            F::one() / x.cosh()
        }
    }

    fn f<F: Float>(x: F) -> F {
        let constant1 = F::from(10).unwrap();
        let constant2 = F::from(100).unwrap();
        let constant3 = F::from(1000).unwrap();

        let constant4 = F::from(0.2).unwrap();
        let constant5 = F::from(0.4).unwrap();
        let constant6 = F::from(0.6).unwrap();

        sech(constant1 * (x - constant4)).powi(2)
            + sech(constant2 * (x - constant5)).powi(4)
            + sech(constant3 * (x - constant6)).powi(6)
    }

    let exact = F::from(0.210_802_736_310_181_77).unwrap();

    Problem {
        id: 21,
        function: f,
        limits: (F::zero(), F::one()),
        exact,
        n: NUM_STEPS,
    }
}

pub fn problem22<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        F::one() / (x.powi(4) + x.powi(2) + F::one())
    }

    let exact = F::from(0.728_102_913_225_581_8).unwrap();

    Problem {
        id: 22,
        function: f,
        limits: (F::zero(), F::one()),
        exact,
        n: NUM_STEPS,
    }
}

pub fn problem23<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        (F::one() / x) * (F::one() / x).sin()
    }

    let constant = F::from(0.1).unwrap();

    let exact = F::from(0.71226).unwrap();

    Problem {
        id: 23,
        function: f,
        limits: (constant, F::one()),
        exact,
        n: NUM_STEPS,
    }
}

pub fn problem24<F: Float + Send + Sum + Sync>() -> Problem<F> {
    fn f<F: Float + Send + Sum + Sync>(x: F) -> F {
        let pi = F::from(PI).unwrap();
        let two = F::one() + F::one();
        let seven = F::from(7).unwrap();

        let n: i32 = 40;

        let sum: F = (1..n)
            .into_par_iter()
            .map(|i| two.powi(-i) * (seven.powi(i) * pi * x).cos())
            .sum();

        sum / pi
    }

    let exact = F::from(-0.0067547455).unwrap();
    let two = F::one() + F::one();

    Problem {
        id: 24,
        function: f,
        limits: (F::zero(), F::one() / two),
        exact,
        n: NUM_STEPS,
    }
}

pub fn problem25<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        let constant = F::from(0.7).unwrap();
        let epsilon = F::from(10e-4).unwrap();

        if x <= constant + epsilon && x >= constant - epsilon {
            F::zero()
        } else {
            (x - constant).abs().ln()
        }
    }

    let exact = F::from(-1.595048).unwrap();

    Problem {
        id: 25,
        function: f,
        limits: (F::zero(), F::one()),
        exact,
        n: NUM_STEPS,
    }
}

pub fn problem26<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        x.cos().exp()
    }
    let two = F::one() + F::one();
    let pi = F::from(PI).unwrap();
    let exact = F::from(7.954_926_521_012_846).unwrap();

    Problem {
        id: 26,
        function: f,
        limits: (F::zero(), two * pi),
        exact,
        n: NUM_STEPS,
    }
}

#[allow(dead_code)]
pub fn problem27<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        let one_over_two = F::one() / (F::one() + F::one());
        let one_over_three = F::one() / (F::one() + F::one() + F::one());

        if x.is_zero() {
            F::zero()
        } else {
            F::one() / (x.powf(one_over_two) + x.powf(one_over_three))
        }
    }
    let exact = F::from(5.0 - 6.0 * 2.0.ln()).unwrap();

    Problem {
        id: 27,
        function: f,
        limits: (F::zero(), F::one()),
        exact,
        n: NUM_STEPS,
    }
}

#[allow(dead_code)]
pub fn problem28<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        let constant = F::from(50).unwrap();

        (-x).exp() * (constant * x).sin()
    }
    let two = F::one() + F::one();
    let pi = F::from(PI).unwrap();
    let exact = F::from(0.019_954_669_277_654_78).unwrap();

    Problem {
        id: 28,
        function: f,
        limits: (F::zero(), two * pi),
        exact,
        n: NUM_STEPS,
    }
}

pub fn problem29<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        let two = F::one() + F::one();
        let threshold = F::from(1.0.exp() - 2.0).unwrap();

        if x > F::zero() && x < threshold {
            F::one() / (x + two)
        } else {
            F::zero()
        }
    }

    let exact = F::from(1.0 - 2.0.ln()).unwrap();

    Problem {
        id: 29,
        function: f,
        limits: (F::zero(), F::one()),
        exact,
        n: NUM_STEPS,
    }
}

pub fn problem30<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        let two = F::one() + F::one();
        let five = F::from(5).unwrap();
        let seven = F::from(7).unwrap();
        let nine = F::from(9).unwrap();

        let constant1 = F::from(1.6).unwrap();
        let constant2 = F::from(4.5).unwrap();

        x.cos() + five * (constant1 * x).cos() - two * (two * x).cos()
            + five * (constant2 * x).cos()
            + seven * (nine * x).cos()
    }

    let exact = F::from(-4.527_569_625_160_672).unwrap();

    let two = F::one() + F::one();
    let seven = F::from(7).unwrap();

    Problem {
        id: 30,
        function: f,
        limits: (two, seven),
        exact,
        n: NUM_STEPS,
    }
}
