//! Gauss-Laguerre quadrature formulas are used to integrate functions $f(x) e^{-x}$ over the positive $x$-axis.
//!
//! With respect to the inner product
//!
//! ```math
//!  \langle f,g \rangle = \int_{0}^{+\infty} f(x) * g(x) * w(x) dx
//! ```
//!
//! the Laguerre polynomials $L_n(x) = e^x \dfrac{\partial^{n} x^n e^{-x}}{\partial x^n}$ for $n > 0$, and $L_0(x) = 1$ form an orthogonal
//! family of polynomials with weight function $w(x) = e^{-x}$ on the positive $x$-axis.
//!
//! The $n$-point Gauss-Laguerre quadrature formula, $GL_n ( f(x) )$, for approximating the integral of $f(x) e^{-x}$ over $\left[0, \infty \right[$,
//! is given by
//! ```math
//! GL_n ( f(x) ) = A_1 f(x_1) + ··· + A_n f(x_n)
//! ```
//! where $xi$ , $i = 1,...,n$, are the zeros of $L_n$ and
//!
//! ```math
//! A_i = \dfrac{n!^2}{ x_i  L_{n-1} (x_i) ^2  }
//! ```
//! for $i = 1,...,n$.

use std::{fmt::Debug, iter::Sum, marker::PhantomData, ops::AddAssign};

use num::{one, Float, One, Zero};
use rayon::iter::{
    IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator, ParallelIterator,
};

use crate::utils::{
    matrix::TridiagonalSymmetricFloatMatrix, orthogonal_polynomials::OrthogonalPolynomial,
};

use super::utils::check_gauss_rule_args;

#[derive(Clone, Debug)]
pub struct Laguerre<F: Float> {
    degree: usize,
    _x: PhantomData<F>,
}

impl<F: Float + Sync + Send + AddAssign + Debug> OrthogonalPolynomial<F> for Laguerre<F> {
    fn new(degree: usize) -> Self {
        Laguerre {
            degree,
            _x: PhantomData,
        }
    }

    fn eval(&self, x: F) -> F {
        if self.degree.is_zero() {
            return F::one();
        }

        if self.degree.is_one() {
            return F::one() - x;
        }

        let mut l_k_1 = F::one(); // L_{k-1}
        let mut l_k = F::one() - x; // L_k

        let mut l = F::nan();

        for k in 2..=self.degree {
            let a = F::from(2 * (k - 1) + 1).unwrap();
            let b = F::from(k - 1).unwrap();
            let c = F::from(k).unwrap();

            l = ((a - x) * l_k - b * l_k_1) / c; // L_{k+1}

            l_k_1 = l_k;
            l_k = l;
        }

        l
    }

    // /// this method used in the process of computing zeros and weights,
    // /// was implemented for some tests, maybe to be removed later
    // /// depending on how OrthogonalPolynomial trait is evolving
    // fn eval_derivative(&self, x: F) -> F {
    //     if x.is_zero() || self.degree.is_zero() {
    //         return zero();
    //     }

    //     let l_n_x = self.eval(x); // L_n(x)
    //     let l_n_1_x = LaguerrePolynomial::new(self.degree - 1).eval(x); // L_{n-1}(x)

    //     let n = F::from(self.degree).unwrap(); // converting n to Float to compute derivative

    //     (n * l_n_x - n * l_n_1_x) / x
    // }

    fn zeros(&self) -> Vec<F> {
        if self.degree.is_zero() {
            return vec![];
        }
        // define the Jacobi matrix (tridiagonal symmetric matrix)

        // we first define the sub-diagonal
        let offdiagonal: Vec<F> = (0..self.degree)
            .into_par_iter()
            .map(|o| F::from(o).unwrap())
            .collect();

        // then the diagonal
        let diagonal: Vec<F> = (0..self.degree)
            .into_par_iter()
            .map(|i| {
                let d = 2 * i + 1;
                F::from(d).unwrap()
            })
            .collect();

        let matrix = TridiagonalSymmetricFloatMatrix::new(diagonal, offdiagonal);

        matrix.eigenvalues()
    }
}

fn roots_laguerre<F: Float + Debug + Sync + Send + AddAssign>(n: usize) -> (Vec<F>, Vec<F>) {
    let l_n: Laguerre<F> = Laguerre::new(n);
    let l_n_plus_1: Laguerre<F> = Laguerre::new(n + 1);

    let zeros = l_n.zeros();

    let n = F::from(n).unwrap();
    let two = F::one() + F::one();

    let weights: Vec<F> = zeros
        .par_iter()
        .map(|x_i| {
            let numerator = *x_i;
            let denominator = (n + one()).powf(two) * l_n_plus_1.eval(*x_i).powf(two);

            numerator / denominator
        })
        .collect();

    let warn = zeros
        .as_slice()
        .into_par_iter()
        .zip(weights.as_slice())
        .any(|(zero, weight)| (*zero).is_nan() || (*weight).is_nan());

    if warn {
        eprintln!(
            "Warning: `n` chosen is too big, some values of Laguerre Polynomials weights or zeros are too small and may underflow!"
        )
    }

    (zeros, weights)
}

pub fn gauss_laguerre_rule<F: Float + Debug + Sync + Send + AddAssign + Sum>(
    f: fn(F) -> F,
    n: usize,
) -> F {
    check_gauss_rule_args(n);
    let (zeros, weights) = roots_laguerre::<F>(n);

    weights
        .into_par_iter()
        .zip(zeros)
        .map(|(w, x)| w * f(x))
        .sum()
}

#[cfg(test)]
mod tests {
    use rayon::iter::IndexedParallelIterator;

    use super::*;
    // use test::Bencher;

    const EPSILON: f64 = 10e-10;

    const L_N_X: &[f64; 17] = &[
        1.0,
        0.0000000000000000E+00,
        -5E-1,
        -6.666_666_666_666_667E-1,
        -6.25E-1,
        -4.666_666_666_666_667E-1,
        -2.569_444_444_444_444E-1,
        -4.047_619_047_619_048E-2,
        1.539_930_555_555_556E-1,
        3.097_442_680_776_014E-1,
        4.189_459_325_396_825E-1,
        4.801_341_790_925_124E-1,
        4.962_122_235_082_305E-1,
        -4.455_729_166_666_667E-1,
        8.5E-1,
        -3.166_666_666_666_667,
        3.433_333_333_333_333E1,
    ];

    // const L_N_X_DERIV: &[f64; 17] = &[
    //     0.0,
    //     -1.0,
    //     -1.0,
    //     -0.5,
    //     0.16666666666666785,
    //     0.7916666666666679,
    //     1.258333333333331,
    //     1.5152777777777793,
    //     1.5557539682539598,
    //     1.4017609126984052,
    //     1.092016644620804,
    //     0.673070712081147,
    //     0.19293653298857372,
    //     -1.1484374999999991,
    //     -0.8749999999999991,
    //     -1.8749999999999991,
    //     11.666666666666785,
    // ];

    const N_VALUES: &[usize; 17] = &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 5, 5, 5, 5];

    const X_VALUES: &[f64; 17] = &[
        1.0E+00, 1.0E+00, 1.0E+00, 1.0E+00, 1.0E+00, 1.0E+00, 1.0E+00, 1.0E+00, 1.0E+00, 1.0E+00,
        1.0E+00, 1.0E+00, 1.0E+00, 0.5E+00, 3.0E+00, 5.0E+00, 1.0E+01,
    ];

    // computed using scipy roots_laguerre
    const FIRST_100_LAGUERRE_ROOTS: [f64; 100] = [
        374.984113,
        355.261312,
        339.435102,
        325.691263,
        313.329534,
        301.985855,
        291.440134,
        281.546328,
        272.20117,
        263.328168,
        254.868629,
        246.776241,
        239.01363,
        231.550068,
        224.359895,
        217.421393,
        210.715973,
        204.22756,
        197.942133,
        191.847369,
        185.93236,
        180.187391,
        174.603761,
        169.17364,
        163.889946,
        158.746249,
        153.736688,
        148.855901,
        144.09897,
        139.461365,
        134.938905,
        130.527723,
        126.224231,
        122.025092,
        117.927199,
        113.927651,
        110.023736,
        106.212911,
        102.492795,
        98.861146,
        95.3158573,
        91.8549433,
        88.4765308,
        85.1788512,
        81.9602322,
        78.8190913,
        75.7539295,
        72.7633254,
        69.8459306,
        67.0004645,
        64.2257101,
        61.5205103,
        58.883764,
        56.314423,
        53.8114889,
        51.3740103,
        49.0010802,
        46.6918335,
        44.4454451,
        42.2611275,
        40.1381291,
        38.0757324,
        36.0732524,
        34.1300351,
        32.245456,
        30.4189188,
        28.6498542,
        26.9377187,
        25.2819939,
        23.6821848,
        22.1378194,
        20.648448,
        19.2136416,
        17.8329919,
        16.5061105,
        15.2326276,
        14.0121922,
        12.8444712,
        11.7291485,
        10.6659251,
        9.65451824,
        8.69466111,
        7.78610238,
        6.92860583,
        6.12195003,
        5.36592799,
        4.66034684,
        4.00502758,
        3.39980483,
        2.84452654,
        2.33905385,
        1.88326083,
        1.47703433,
        1.12027384,
        0.812891284,
        0.554810938,
        0.345969181,
        0.186314102,
        0.075803612,
        0.014386147,
    ];

    const FIRST_100_LAGUERRE_WEIGHTS: [f64; 100] = [
        3.24656516e-162,
        8.90503141e-154,
        5.6260373e-147,
        4.64686301e-141,
        9.8824946e-136,
        7.71361149e-131,
        2.73996547e-126,
        5.11064048e-122,
        5.53964175e-118,
        3.7620973e-114,
        1.69596926e-110,
        5.31213273e-107,
        1.19948442e-103,
        2.01265479e-100,
        2.57396151e-97,
        2.56339222e-94,
        2.02484835e-91,
        1.28896375e-88,
        6.70480123e-86,
        2.88487552e-83,
        1.03790059e-80,
        3.15247255e-78,
        8.15375261e-76,
        1.80985953e-73,
        3.47185478e-71,
        5.79263076e-69,
        8.45495649e-67,
        1.08536745e-64,
        1.23137657e-62,
        1.24023504e-60,
        1.11356697e-58,
        8.94733725e-57,
        6.4562691e-55,
        4.19774702e-53,
        2.46681069e-51,
        1.3139812e-49,
        6.36127057e-48,
        2.8060375e-46,
        1.1304819e-44,
        4.168864e-43,
        1.41013923e-41,
        4.38379098e-40,
        1.254833e-38,
        3.3130638e-37,
        8.08163234e-36,
        1.8242002e-34,
        3.81586014e-33,
        7.40742896e-32,
        1.33620942e-30,
        2.24264627e-29,
        3.50627215e-28,
        5.11236982e-27,
        6.95919878e-26,
        8.85323177e-25,
        1.05359406e-23,
        1.17402171e-22,
        1.22600701e-21,
        1.20084697e-20,
        1.10409452e-19,
        9.53624945e-19,
        7.7430891e-18,
        5.91443247e-17,
        4.2526129e-16,
        2.88011506e-15,
        1.83835018e-14,
        1.10649802e-13,
        6.28352496e-13,
        3.36821417e-12,
        1.70506256e-11,
        8.15479892e-11,
        3.6863292e-10,
        1.57560032e-09,
        6.36971137e-09,
        2.43642586e-08,
        8.8200584e-08,
        3.02263874e-07,
        9.8083359e-07,
        3.01426749e-06,
        8.77430976e-06,
        2.41957523e-05,
        6.32108705e-05,
        0.000156452074,
        0.000366854837,
        0.000814871592,
        0.00171431974,
        0.00341497999,
        0.006438951,
        0.0114854424,
        0.0193678281,
        0.0308463086,
        0.0463401336,
        0.0655510093,
        0.0870966385,
        0.108314112,
        0.125407091,
        0.13404334,
        0.130356613,
        0.112115103,
        0.0796767462,
        0.0363926059,
    ];

    #[test]
    fn test_laguerre_polynomial_zeros() {
        const EPSILON: f64 = 10e-5;

        let n = 100;
        let lag: Laguerre<f64> = Laguerre::new(n);

        let lag_zeros = lag.zeros();

        FIRST_100_LAGUERRE_ROOTS
            .into_par_iter()
            .zip(lag_zeros)
            .for_each(|(test_zero, zero)| assert!((test_zero - zero).abs() < EPSILON))
    }

    #[test]
    fn test_laguerre_polynomial_weights() {
        const EPSILON: f64 = 10e-5;

        let n = 100;
        let (_, weights) = roots_laguerre::<f64>(n);

        FIRST_100_LAGUERRE_WEIGHTS
            .into_par_iter()
            .zip(weights)
            .for_each(|(test_weight, weight)| assert!((test_weight - weight).abs() < EPSILON))
    }

    #[test]
    fn test_eval_laguerre() {
        for ((&ln_test, &n), &x) in L_N_X.iter().zip(N_VALUES).zip(X_VALUES) {
            let lag: Laguerre<f64> = Laguerre::new(n);

            let ln = lag.eval(x);

            let is_close = (ln - ln_test).abs() < EPSILON;

            assert!(is_close)
        }
    }

    // #[test]
    // fn test_eval_laguerre_derivative() {
    //     for ((&dln_test, &n), &x) in L_N_X_DERIV.iter().zip(N_VALUES).zip(X_VALUES) {
    //         let lag: LaguerrePolynomial<f64> = LaguerrePolynomial::new(n);

    //         let dln = lag.eval_derivative(x);

    //         let is_close = (dln - dln_test).abs() < EPSILON;

    //         assert!(is_close)
    //     }
    // }

    // #[bench]
    // fn bench_roots_laguerre(bencher: &mut Bencher) {
    //     let n: usize = 1_000;

    //     bencher.iter(|| {
    //         roots_laguerre::<f64>(n);
    //     })
    // }
}
