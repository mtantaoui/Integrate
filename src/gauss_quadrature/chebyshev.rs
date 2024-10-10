use std::iter::Sum;
use std::{f64::consts::PI, marker::PhantomData};

use std::fmt::Debug;
use std::ops::AddAssign;

use num::{one, Float, Zero};
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};

use crate::utils::orthogonal_polynomials::OrthogonalPolynomial;

use super::utils::check_gauss_rule_args;

#[derive(Clone, Debug)]
struct ChebyshevFirstKind<F: Float> {
    degree: usize,
    _x: PhantomData<F>,
}

#[derive(Clone, Debug)]
struct ChebyshevSecondKind<F: Float> {
    degree: usize,
    _x: PhantomData<F>,
}

fn roots_first_kind_chebyshev<F: Float + Debug + Sync + Send + AddAssign>(
    n: usize,
) -> (Vec<F>, Vec<F>) {
    let t_n: ChebyshevFirstKind<F> = ChebyshevFirstKind::new(n);
    let zeros = t_n.zeros();

    let pi = F::from(PI).unwrap();
    let n = F::from(n).unwrap();

    let weights = vec![pi / n; t_n.degree];

    let warn = zeros
        .as_slice()
        .into_par_iter()
        .zip(weights.as_slice())
        .any(|(zero, weight)| (*zero).is_nan() || (*weight).is_nan());

    if warn {
        eprintln!(
        "Warning: `n` chosen is too big, some values of Chebyshev First Kind Polynomials weights or zeros are too small and may underflow!"
    );
    }

    (zeros, weights)
}

fn roots_second_kind_chebyshev<F: Float + Debug + Sync + Send + AddAssign>(
    n: usize,
) -> (Vec<F>, Vec<F>) {
    let u_n: ChebyshevSecondKind<F> = ChebyshevSecondKind::new(n);
    let zeros = u_n.zeros();

    let n = F::from(u_n.degree).unwrap();
    let pi = F::from(PI).unwrap();

    let weights: Vec<F> = (1..=u_n.degree)
        .into_par_iter()
        .map(|i| {
            let i = F::from(i).unwrap();

            let numer = i * pi;
            let denom = n + one();

            let angle = numer / denom;

            let term1 = angle.sin().powi(2);
            let term2 = pi / (n + one());

            term1 * term2
        })
        .collect();

    let warn = zeros
        .as_slice()
        .into_par_iter()
        .zip(weights.as_slice())
        .any(|(zero, weight)| (*zero).is_nan() || (*weight).is_nan());

    if warn {
        eprintln!(
            "Warning: `n` chosen is too big, some values of Chebyshev Second Kind Polynomials weights or zeros are too small and may underflow!"
        );
    }

    (zeros, weights)
}

pub fn gauss_first_kind_chebyshev_rule<F: Float + Debug + Sync + Send + AddAssign + Sum>(
    f: fn(F) -> F,
    n: usize,
) -> F {
    check_gauss_rule_args(n);

    let (zeros, weights) = roots_first_kind_chebyshev::<F>(n);

    weights
        .into_par_iter()
        .zip(zeros)
        .map(|(w, x)| w * f(x))
        .sum()
}

pub fn gauss_second_kind_chebyshev_rule<F: Float + Debug + Sync + Send + AddAssign + Sum>(
    f: fn(F) -> F,
    n: usize,
) -> F {
    check_gauss_rule_args(n);

    let (zeros, weights) = roots_second_kind_chebyshev::<F>(n);

    weights
        .into_par_iter()
        .zip(zeros)
        .map(|(w, x)| w * f(x))
        .sum()
}

impl<F: Float + Debug + AddAssign + Send + Sync> OrthogonalPolynomial<F> for ChebyshevFirstKind<F> {
    fn new(degree: usize) -> Self {
        ChebyshevFirstKind {
            degree,
            _x: PhantomData,
        }
    }

    fn eval(&self, x: F) -> F {
        let theta = x.acos();
        let n = F::from(self.degree).unwrap();

        (n * theta).cos()
    }

    fn zeros(&self) -> Vec<F> {
        if self.degree.is_zero() {
            return vec![];
        }

        let n = F::from(self.degree).unwrap();
        let pi = F::from(PI).unwrap();
        let two = F::one() + F::one();

        let zeros: Vec<F> = (1..=self.degree)
            .into_par_iter()
            .map(|i| {
                let i = F::from(i).unwrap();

                let numer = (two * i - one()) * pi;
                let denom = two * n;

                let angle = numer / denom;

                angle.cos()
            })
            .collect();

        zeros
    }
}

impl<F: Float + Debug + AddAssign + Send + Sync> OrthogonalPolynomial<F>
    for ChebyshevSecondKind<F>
{
    fn new(degree: usize) -> Self {
        ChebyshevSecondKind {
            degree,
            _x: PhantomData,
        }
    }

    fn eval(&self, x: F) -> F {
        let theta = x.acos();
        let n = F::from(self.degree).unwrap();

        let numer = ((n + one()) * theta).sin();
        let denom = theta.sin();

        numer / denom
    }

    fn zeros(&self) -> Vec<F> {
        if self.degree.is_zero() {
            return vec![];
        }

        let n = F::from(self.degree).unwrap();
        let pi = F::from(PI).unwrap();

        let zeros: Vec<F> = (1..=self.degree)
            .into_par_iter()
            .map(|i| {
                let i = F::from(i).unwrap();

                let numer = i * pi;
                let denom = n + one();

                let angle = numer / denom;

                angle.cos()
            })
            .collect();

        zeros
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, FRAC_PI_8, PI};

    use crate::{
        gauss_quadrature::chebyshev::{
            gauss_first_kind_chebyshev_rule, gauss_second_kind_chebyshev_rule,
            roots_first_kind_chebyshev, roots_second_kind_chebyshev, ChebyshevFirstKind,
            ChebyshevSecondKind,
        },
        utils::orthogonal_polynomials::OrthogonalPolynomial,
    };

    const EPSILON: f64 = 10e-6;

    const T1_ZEROS: [f64; 1] = [6.123_233_995_736_766E-17];
    const T2_ZEROS: [f64; 2] = [-0.7071067811865475, 0.7071067811865475];
    const T4_ZEROS: [f64; 4] = [
        -0.9238795325112867,
        -0.3826834323650897,
        0.3826834323650898,
        0.9238795325112867,
    ];
    const T8_ZEROS: [f64; 8] = [
        -0.9807852804032304,
        -0.8314696123025453,
        -0.555_570_233_019_602,
        -0.1950903220161282,
        0.1950903220161283,
        0.5555702330196023,
        0.8314696123025452,
        0.9807852804032304,
    ];
    const T16_ZEROS: [f64; 16] = [
        -0.9951847266721968,
        -0.9569403357322088,
        -0.8819212643483549,
        -0.773_010_453_362_737,
        -0.6343932841636454,
        -0.4713967368259977,
        -0.2902846772544622,
        -9.801_714_032_956_065E-2,
        9.801_714_032_956_077E-2,
        0.2902846772544623,
        0.4713967368259978,
        0.6343932841636455,
        0.773_010_453_362_737,
        0.881_921_264_348_355,
        0.9569403357322088,
        0.9951847266721969,
    ];

    const T1_WEIGHTS: [f64; 1] = [PI];
    const T2_WEIGHTS: [f64; 2] = [FRAC_PI_2, FRAC_PI_2];
    const T4_WEIGHTS: [f64; 4] = [FRAC_PI_4, FRAC_PI_4, FRAC_PI_4, FRAC_PI_4];
    const T8_WEIGHTS: [f64; 8] = [
        FRAC_PI_8, FRAC_PI_8, FRAC_PI_8, FRAC_PI_8, FRAC_PI_8, FRAC_PI_8, FRAC_PI_8, FRAC_PI_8,
    ];
    const T16_WEIGHTS: [f64; 16] = [
        0.1963495408493621,
        0.1963495408493621,
        0.1963495408493621,
        0.1963495408493621,
        0.1963495408493621,
        0.1963495408493621,
        0.1963495408493621,
        0.1963495408493621,
        0.1963495408493621,
        0.1963495408493621,
        0.1963495408493621,
        0.1963495408493621,
        0.1963495408493621,
        0.1963495408493621,
        0.1963495408493621,
        0.1963495408493621,
    ];

    const U1_ZEROS: [f64; 1] = [6.123_233_995_736_766E-17];
    const U2_ZEROS: [f64; 2] = [-0.4999999999999998, 0.5000000000000001];
    const U4_ZEROS: [f64; 4] = [
        -0.8090169943749473,
        -0.3090169943749473,
        0.3090169943749475,
        0.8090169943749475,
    ];
    const U8_ZEROS: [f64; 8] = [
        -0.9396926207859083,
        -0.7660444431189779,
        -0.4999999999999998,
        -0.1736481776669303,
        0.1736481776669304,
        0.5000000000000001,
        0.766_044_443_118_978,
        0.9396926207859084,
    ];
    const U16_ZEROS: [f64; 16] = [
        -0.9829730996839018,
        -0.9324722294043557,
        -0.850_217_135_729_614,
        -0.7390089172206593,
        -0.6026346363792563,
        -0.4457383557765379,
        -0.2736629900720829,
        -9.226_835_946_330_189E-2,
        9.226_835_946_330_202E-2,
        0.273_662_990_072_083,
        0.4457383557765384,
        0.6026346363792564,
        0.7390089172206592,
        0.8502171357296141,
        0.9324722294043558,
        0.9829730996839018,
    ];

    const U1_WEIGHTS: [f64; 1] = [FRAC_PI_2];

    const U2_WEIGHTS: [f64; 2] = [0.7853981633974484, 0.7853981633974481];

    const U4_WEIGHTS: [f64; 4] = [
        0.2170787134227061,
        0.5683194499747424,
        0.5683194499747423,
        0.217_078_713_422_706,
    ];

    const U8_WEIGHTS: [f64; 8] = [
        4.083_294_770_910_712E-2,
        0.1442256007956728,
        0.2617993877991495,
        0.338_540_227_093_519,
        0.338_540_227_093_519,
        0.2617993877991494,
        0.1442256007956728,
        4.083_294_770_910_708E-2,
    ];
    const U16_WEIGHTS: [f64; 16] = [
        6.239_551_412_252_139E-3,
        2.411_551_965_623_602E-2,
        5.121_365_616_644_202E-2,
        8.387_420_745_120_884E-2,
        0.1176861850811667,
        0.1480830941187536,
        0.170_959_663_563_356,
        0.1832262859480331,
        0.1832262859480331,
        0.170_959_663_563_356,
        0.1480830941187536,
        0.1176861850811666,
        8.387_420_745_120_885E-2,
        5.121_365_616_644_197E-2,
        2.411_551_965_623_597E-2,
        6.239_551_412_252_137E-3,
    ];

    const DEGREE: [usize; 13] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];

    const EVAL1: [f64; 13] = [
        1.0, 0.8, 0.28, -0.352, -0.8432, -0.99712, -0.752192, -0.206387, 0.421972, 0.881543,
        0.988497, 0.700051, 0.131586,
    ];
    const EVAL2: [f64; 13] = [
        1.0, 1.6, 1.56, 0.896, -0.1264, -1.09824, -1.63078, -1.51101, -0.786839, 0.252072, 1.19015,
        1.65217, 1.45333,
    ];

    #[test]
    fn test_chebyshev_first_kind_zeros() {
        let t1: ChebyshevFirstKind<f64> = ChebyshevFirstKind::new(1);
        let t2: ChebyshevFirstKind<f64> = ChebyshevFirstKind::new(2);
        let t4: ChebyshevFirstKind<f64> = ChebyshevFirstKind::new(4);
        let t8: ChebyshevFirstKind<f64> = ChebyshevFirstKind::new(8);
        let t16: ChebyshevFirstKind<f64> = ChebyshevFirstKind::new(16);

        let mut t1_zeros = t1.zeros();
        t1_zeros.reverse();

        let mut t2_zeros = t2.zeros();
        t2_zeros.reverse();

        let mut t4_zeros = t4.zeros();
        t4_zeros.reverse();

        let mut t8_zeros = t8.zeros();
        t8_zeros.reverse();

        let mut t16_zeros = t16.zeros();
        t16_zeros.reverse();

        let t1_test = t1_zeros
            .iter()
            .zip(T1_ZEROS)
            .all(|(value, test_value)| (value - test_value).abs() <= EPSILON);
        assert!(t1_test);

        let t2_test = t2_zeros
            .iter()
            .zip(T2_ZEROS)
            .all(|(value, test_value)| (value - test_value).abs() <= EPSILON);
        assert!(t2_test);

        let t4_test = t4_zeros
            .iter()
            .zip(T4_ZEROS)
            .all(|(value, test_value)| (value - test_value).abs() <= EPSILON);
        assert!(t4_test);

        let t8_test = t8_zeros
            .iter()
            .zip(T8_ZEROS)
            .all(|(value, test_value)| (value - test_value).abs() <= EPSILON);
        assert!(t8_test);

        let t16_test = t16_zeros
            .iter()
            .zip(T16_ZEROS)
            .all(|(value, test_value)| (value - test_value).abs() <= EPSILON);
        assert!(t16_test)
    }

    #[test]
    fn test_chebyshev_first_kind_weights() {
        let (_, t1_weights) = roots_first_kind_chebyshev::<f64>(1);
        let (_, t2_weights) = roots_first_kind_chebyshev::<f64>(2);
        let (_, t4_weights) = roots_first_kind_chebyshev::<f64>(4);
        let (_, t8_weights) = roots_first_kind_chebyshev::<f64>(8);
        let (_, t16_weights) = roots_first_kind_chebyshev::<f64>(16);

        let t1_test = t1_weights
            .iter()
            .zip(T1_WEIGHTS)
            .all(|(value, test_value)| (value - test_value).abs() <= EPSILON);
        assert!(t1_test);

        let t2_test = t2_weights
            .iter()
            .zip(T2_WEIGHTS)
            .all(|(value, test_value)| (value - test_value).abs() <= EPSILON);
        assert!(t2_test);

        let t4_test = t4_weights
            .iter()
            .zip(T4_WEIGHTS)
            .all(|(value, test_value)| (value - test_value).abs() <= EPSILON);
        assert!(t4_test);

        let t8_test = t8_weights
            .iter()
            .zip(T8_WEIGHTS)
            .all(|(value, test_value)| (value - test_value).abs() <= EPSILON);
        assert!(t8_test);

        let t16_test = t16_weights
            .iter()
            .zip(T16_WEIGHTS)
            .all(|(value, test_value)| (value - test_value).abs() <= EPSILON);
        assert!(t16_test)
    }

    #[test]
    fn test_eval_chebyshev_first_kind_zeros() {
        for (i, test_value) in DEGREE.into_iter().zip(EVAL1) {
            let t_n: ChebyshevFirstKind<f64> = ChebyshevFirstKind::new(i);

            let computed = t_n.eval(0.8);

            assert!((computed - test_value).abs() < EPSILON)
        }
    }

    #[test]
    fn test_chebyshev_second_kind_zeros() {
        let u1: ChebyshevSecondKind<f64> = ChebyshevSecondKind::new(1);
        let u2: ChebyshevSecondKind<f64> = ChebyshevSecondKind::new(2);
        let u4: ChebyshevSecondKind<f64> = ChebyshevSecondKind::new(4);
        let u8: ChebyshevSecondKind<f64> = ChebyshevSecondKind::new(8);
        let u16: ChebyshevSecondKind<f64> = ChebyshevSecondKind::new(16);

        let mut u1_zeros = u1.zeros();
        u1_zeros.reverse();

        let mut u2_zeros = u2.zeros();
        u2_zeros.reverse();

        let mut u4_zeros = u4.zeros();
        u4_zeros.reverse();

        let mut u8_zeros = u8.zeros();
        u8_zeros.reverse();

        let mut u16_zeros = u16.zeros();
        u16_zeros.reverse();

        let u1_test = u1_zeros
            .iter()
            .zip(U1_ZEROS)
            .all(|(value, test_value)| (value - test_value).abs() <= EPSILON);
        assert!(u1_test);

        let u2_test = u2_zeros
            .iter()
            .zip(U2_ZEROS)
            .all(|(value, test_value)| (value - test_value).abs() <= EPSILON);
        assert!(u2_test);

        let u4_test = u4_zeros
            .iter()
            .zip(U4_ZEROS)
            .all(|(value, test_value)| (value - test_value).abs() <= EPSILON);
        assert!(u4_test);

        let u8_test = u8_zeros
            .iter()
            .zip(U8_ZEROS)
            .all(|(value, test_value)| (value - test_value).abs() <= EPSILON);
        assert!(u8_test);

        let u16_test = u16_zeros
            .iter()
            .zip(U16_ZEROS)
            .all(|(value, test_value)| (value - test_value).abs() <= EPSILON);
        assert!(u16_test)
    }

    #[test]
    fn test_chebyshev_second_kind_weights() {
        let (_, u1_weights) = roots_second_kind_chebyshev::<f64>(1);
        let (_, u2_weights) = roots_second_kind_chebyshev::<f64>(2);
        let (_, u4_weights) = roots_second_kind_chebyshev::<f64>(4);
        let (_, u8_weights) = roots_second_kind_chebyshev::<f64>(8);
        let (_, u16_weights) = roots_second_kind_chebyshev::<f64>(16);

        let u1_test = u1_weights
            .iter()
            .zip(U1_WEIGHTS)
            .all(|(value, test_value)| (value - test_value).abs() <= EPSILON);
        assert!(u1_test);

        let u2_test = u2_weights
            .iter()
            .zip(U2_WEIGHTS)
            .all(|(value, test_value)| (value - test_value).abs() <= EPSILON);
        assert!(u2_test);

        let u4_test = u4_weights
            .iter()
            .zip(U4_WEIGHTS)
            .all(|(value, test_value)| (value - test_value).abs() <= EPSILON);
        assert!(u4_test);

        let u8_test = u8_weights
            .iter()
            .zip(U8_WEIGHTS)
            .all(|(value, test_value)| (value - test_value).abs() <= EPSILON);
        assert!(u8_test);

        let u16_test = u16_weights
            .iter()
            .zip(U16_WEIGHTS)
            .all(|(value, test_value)| (value - test_value).abs() <= EPSILON);
        assert!(u16_test)
    }

    #[test]
    fn test_eval_chebyshev_second_kind_zeros() {
        for (i, test_value) in DEGREE.into_iter().zip(EVAL2) {
            let u_n: ChebyshevSecondKind<f64> = ChebyshevSecondKind::new(i);

            let computed = u_n.eval(0.8);

            assert!((computed - test_value).abs() < EPSILON)
        }
    }

    // Test the numerical integration of cos(1000 x) over the range [-1,1]
    // for varying number of Gauss-Chebyshev First Kind quadrature nodes l.
    // exact value of the numerical integration is 0.002 * sin(1000)
    // The fact that only twelve digits of accuracy are obtained is due to the
    // condition number of the summation.
    #[test]
    fn test_chebyshev_first_kind_rule() {
        let exact: f64 = 0.002 * (1000.0_f64).sin();

        fn f(x: f64) -> f64 {
            (1000.0 * x).cos() * (1.0 - x.powi(2)).sqrt()
        }

        println!("Integral Exact Value: {}", exact);

        for l in (540..=700_usize).step_by(20) {
            // Gauss-Legendre rule using glpair function
            let integral: f64 = gauss_first_kind_chebyshev_rule(f, l);

            println!(
                "number of nodes: {} \t Gauss-Chebyshev First Kind Integral: {}",
                l, integral
            );

            assert!(integral - exact < EPSILON);
        }
    }

    // Test the numerical integration of cos(1000 x) over the range [-1,1]
    // for varying number of Gauss-Chebyshev Second Kind quadrature nodes l.
    // exact value of the numerical integration is 0.002 * sin(1000)
    // The fact that only twelve digits of accuracy are obtained is due to the
    // condition number of the summation.
    #[test]
    fn test_chebyshev_second_kind_rule() {
        let exact: f64 = 0.002 * (1000.0_f64).sin();

        fn f(x: f64) -> f64 {
            (1000.0 * x).cos() / (1.0 - x.powi(2)).sqrt()
        }

        println!("Integral Exact Value: {}", exact);

        for l in (540..=700_usize).step_by(20) {
            // Gauss-Legendre rule using glpair function
            let integral: f64 = gauss_second_kind_chebyshev_rule(f, l);

            println!(
                "number of nodes: {} \t Gauss-Chebyshev Second Kind Integral: {}",
                l, integral
            );

            assert!((integral - exact).abs() < EPSILON);
        }
    }
}
