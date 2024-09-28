use std::{f64::consts::PI, i32, iter::Sum};

use num::Float;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

const EPSILON: f64 = 10e-6;
const NUM_STEPS: usize = 100_000;

#[derive(Debug, Clone)]
pub enum Methods {
    Rectangle,
    Trapezoidal,
    Newton3Over8,
    Simpson,
}

impl Methods {
    pub fn iter() -> core::array::IntoIter<Methods, 4> {
        [
            Methods::Rectangle,
            Methods::Trapezoidal,
            Methods::Newton3Over8,
            Methods::Simpson,
        ]
        .into_iter()
    }
}

#[derive(Clone)]
pub struct Problem<F: Float> {
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

fn problem01<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        x.exp()
    }

    Problem {
        function: f,
        limits: (F::zero(), F::one()),
        exact: F::one().exp() - F::one(),
        n: NUM_STEPS,
    }
}

fn problem02<F: Float>() -> Problem<F> {
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
        function: f,
        limits: (F::zero(), F::one()),
        exact: exact,
        n: NUM_STEPS,
    }
}

fn problem03<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        x.sqrt()
    }

    let exact = F::from(2.0 / 3.0).unwrap();

    Problem {
        function: f,
        limits: (F::zero(), F::one()),
        exact: exact,
        n: NUM_STEPS,
    }
}

fn problem04<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        let constant = F::from(0.92).unwrap();
        constant * x.cosh() - x.cos()
    }

    let exact = F::from(1.84 * 1.0.sinh() - 2.0 * 1.0.sin()).unwrap();

    Problem {
        function: f,
        limits: (-F::one(), F::one()),
        exact: exact,
        n: NUM_STEPS,
    }
}

fn problem05<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        let constant = F::from(0.9).unwrap();
        F::one() / (x.powi(4) + x.powi(2) + constant)
    }

    let exact = F::from(1.5822329637296729331).unwrap();

    Problem {
        function: f,
        limits: (-F::one(), F::one()),
        exact: exact,
        n: NUM_STEPS,
    }
}

fn problem06<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        let constant = F::from(0.5).unwrap();
        (x + constant).abs().sqrt()
    }

    let exact = (2.0.sqrt() + 3.0 * 6.0.sqrt()) / 6.0;
    let exact = F::from(exact).unwrap();

    Problem {
        function: f,
        limits: (-F::one(), F::one()),
        exact: exact,
        n: NUM_STEPS,
    }
}

fn problem07<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        F::one() / x.sqrt()
    }

    let constant = F::from(10_000).unwrap();
    let two = F::one() + F::one();

    let exact = F::one() + F::one() - two * (F::one() / constant).sqrt();

    Problem {
        function: f,
        limits: (F::one() / constant, F::one()),
        exact,
        n: NUM_STEPS,
    }
}

fn problem08<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        F::one() / (F::one() + x.powi(4))
    }

    let exact = F::from(0.86697298733991103757).unwrap();

    Problem {
        function: f,
        limits: (F::zero(), F::one()),
        exact,
        n: NUM_STEPS,
    }
}

fn problem09<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        let two = F::one() + F::one();
        let ten = F::from(10).unwrap();
        let pi = F::from(PI).unwrap();

        two / (two + (ten * pi * x).sin())
    }

    let exact = F::from(2.0 / 3.0.sqrt()).unwrap();

    Problem {
        function: f,
        limits: (F::zero(), F::one()),
        exact,
        n: NUM_STEPS,
    }
}

fn problem10<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        F::one() / (F::one() + x)
    }

    let exact = F::from(2.0.ln()).unwrap();

    Problem {
        function: f,
        limits: (F::zero(), F::one()),
        exact,
        n: NUM_STEPS,
    }
}

fn problem11<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        F::one() / (F::one() + x.exp())
    }

    let exact = F::from(((2.0 * 1.0.exp()) / (1.0 + 1.0.exp())).ln()).unwrap();

    Problem {
        function: f,
        limits: (F::zero(), F::one()),
        exact,
        n: NUM_STEPS,
    }
}

fn problem12<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        x / (x.exp() - F::one())
    }

    let exact = F::from(0.77750463411224827642).unwrap();

    Problem {
        function: f,
        limits: (F::epsilon(), F::one()),
        exact,
        n: NUM_STEPS,
    }
}

fn problem13<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        x.sin() / x
    }

    let exact = F::from(1.6583475942188740493).unwrap();

    let ten = F::from(10).unwrap();

    Problem {
        function: f,
        limits: (F::epsilon(), ten),
        exact,
        n: NUM_STEPS,
    }
}

fn problem14<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        let constant = F::from(50).unwrap();
        let pi = F::from(PI).unwrap();
        constant.sqrt() * (-constant * pi * x * x).exp()
    }

    let exact = F::from(0.500000211166).unwrap();

    let ten = F::from(10).unwrap();

    Problem {
        function: f,
        limits: (F::zero(), ten),
        exact,
        n: NUM_STEPS,
    }
}

fn problem15<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        let constant = F::from(25).unwrap();
        constant * (-constant * x).exp()
    }

    let exact = F::from(1.0 - (-250.0).exp()).unwrap();

    let ten = F::from(10).unwrap();

    Problem {
        function: f,
        limits: (F::zero(), ten),
        exact,
        n: NUM_STEPS,
    }
}

fn problem16<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        let constant1 = F::from(50).unwrap();
        let constant2 = F::from(2500).unwrap();

        let pi = F::from(PI).unwrap();

        constant1 / (pi * (constant2 * x * x + F::one()))
    }

    let exact = F::from(0.49936338107645674464).unwrap();

    let ten = F::from(10).unwrap();

    Problem {
        function: f,
        limits: (F::zero(), ten),
        exact,
        n: NUM_STEPS,
    }
}

fn problem17<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        let constant = F::from(50).unwrap();

        let pi = F::from(PI).unwrap();

        (constant * pi * x).sin().powi(2)
    }

    let exact = F::from(0.5).unwrap();

    Problem {
        function: f,
        limits: (F::zero(), F::one()),
        exact,
        n: NUM_STEPS,
    }
}

fn problem18<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        x / (x.exp() + F::one())
    }

    let exact = F::from(0.17055734950243820437).unwrap();

    Problem {
        function: f,
        limits: (F::zero(), F::one()),
        exact,
        n: NUM_STEPS,
    }
}

fn problem19<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        x.ln()
    }

    fn anti_derivative_f<F: Float>(x: F) -> F {
        x * x.ln() - x
    }

    let constant = F::one() / F::from(10_000).unwrap();
    let exact = anti_derivative_f(F::one()) - anti_derivative_f(constant);

    Problem {
        function: f,
        limits: (constant, F::one()),
        exact,
        n: NUM_STEPS,
    }
}

fn problem20<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        let constant = F::from(1.005).unwrap();
        F::one() / (x.powi(2) + constant)
    }

    let exact = F::from(1.5643964440690497731).unwrap();

    Problem {
        function: f,
        limits: (-F::one(), F::one()),
        exact,
        n: NUM_STEPS,
    }
}

fn problem21<F: Float>() -> Problem<F> {
    // evaluates the hyperbolic secant, while avoiding COSH overflow.
    fn sech<F: Float>(x: F) -> F {
        let log_huge = F::from(80.0).unwrap();

        let value = if log_huge < x.abs() {
            F::zero()
        } else {
            F::one() / x.cosh()
        };

        return value;
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

    let exact = F::from(0.21080273631018169851).unwrap();

    Problem {
        function: f,
        limits: (F::zero(), F::one()),
        exact,
        n: NUM_STEPS,
    }
}

fn problem22<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        F::one() / (x.powi(4) + x.powi(2) + F::one())
    }

    let exact = F::from(0.72810291322558188550).unwrap();

    Problem {
        function: f,
        limits: (F::zero(), F::one()),
        exact,
        n: NUM_STEPS,
    }
}

fn problem23<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        (F::one() / x) * (F::one() / x).sin()
    }

    let constant = F::from(0.1).unwrap();

    let exact = F::from(0.71226).unwrap();

    Problem {
        function: f,
        limits: (constant, F::one()),
        exact,
        n: NUM_STEPS,
    }
}

fn problem24<F: Float + Send + Sum + Sync>() -> Problem<F> {
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
        function: f,
        limits: (F::zero(), F::one() / two),
        exact,
        n: NUM_STEPS,
    }
}

fn problem25<F: Float>() -> Problem<F> {
    fn f<F: Float>(x: F) -> F {
        let constant = F::from(0.7).unwrap();
        let epsilon = F::from(10e-6).unwrap();

        if x < constant + epsilon && x > constant - epsilon {
            F::zero()
        } else {
            (x - constant).abs().ln()
        }
    }

    let exact = F::from(-1.6108).unwrap();

    Problem {
        function: f,
        limits: (F::zero(), F::one()),
        exact,
        n: NUM_STEPS,
    }
}

pub fn problems_vec<F: Float + Send + Sum + Sync>() -> Vec<Problem<F>> {
    vec![
        problem01(),
        problem02(),
        problem03(),
        problem04(),
        problem05(),
        problem06(),
        problem07(),
        problem08(),
        problem09(),
        problem10(),
        problem11(),
        problem12(),
        problem13(),
        problem14(),
        problem15(),
        problem16(),
        problem17(),
        problem18(),
        problem19(),
        problem20(),
        problem21(),
        problem22(),
        problem23(),
        problem24(),
        // problem25(),
    ]
}
