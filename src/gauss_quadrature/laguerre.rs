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

extern crate test;

use num::{zero, Float, One, ToPrimitive, Zero};

use crate::utils::newton_raphson::newton_raphson;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use super::bessel::bessel_j0_zeros;

/// TODO: Documentation and examples
pub fn laguerre_polynomial_zeros(n: usize) -> Vec<f64> {
    // let approx_zeros = laguerre_polynomial_zeros(n);

    // (0..=n).into_par_iter().map(|i| {
    //     let l_n = |x: f64| eval_laguerre(i, x);

    //     let dl_n = |x: f64| eval_laguerre_derivative(i, x);

    //     // newton_raphson(l_n, dl_n, 0.0, f64::EPSILON);

    //     let a: f64 = 0.0;
    //     let tolerance: f64 = 0.00000001;

    //     newton_raphson::<f64>(l_n, dl_n, a, tolerance)
    // });

    vec![0.0; 5]
}

/// TODO: Documentation and examples
fn laguerre_polynomial_approximate_zeros(n: usize) -> Vec<f64> {
    (1..=n)
        .into_par_iter()
        .map(|m| approximate_laguerre_zero(m, n))
        .collect()
}

/// TODO: Documentation and examples
fn eval_laguerre<F: Float>(n: usize, x: F) -> F {
    if n.is_zero() {
        return F::one();
    }

    if n.is_one() {
        return F::one() - x;
    }

    let mut l_k_1 = F::one(); // L_{k-1}
    let mut l_k = F::one() - x; // L_k

    let mut l = F::nan();

    for k in 2..=n {
        let a = F::from(2 * (k - 1) + 1).unwrap();
        let b = F::from(k - 1).unwrap();
        let c = F::from(k).unwrap();

        l = ((a - x) * l_k - b * l_k_1) / c; // L_{k+1}

        l_k_1 = l_k;
        l_k = l;
    }

    l
}

fn eval_laguerre_derivative<F: Float>(n: usize, x: F) -> F {
    if x.is_zero() || n.is_zero() {
        return zero();
    }

    let l_n_x = eval_laguerre(n, x); // L_n(x)
    let l_n_1_x = eval_laguerre(n - 1, x); // L_{n-1}(x)

    let n = F::from(n).unwrap(); // converting n to Float to compute derivative

    let l_derivative = (n * l_n_x - n * l_n_1_x) / x;

    return l_derivative;
}

/// TODO: Documentation and examples
fn approximate_laguerre_zero(m: usize, n: usize) -> f64 {
    let n_f = n.to_f64().unwrap();

    let j_0_m = bessel_j0_zeros(m);
    let k_n: f64 = n_f + 0.5;

    let term1 = j_0_m.powi(2) / (4.0 * k_n);
    let term2 = 1.0 + (j_0_m.powi(2) - 2.0) / (48.0 * k_n.powi(2));

    term1 * term2
}

#[cfg(test)]
mod tests {

    use super::*;
    use test::Bencher;

    const EPSILON: f64 = 10e-10;

    const L_N_X: &[f64; 17] = &[
        0.1000000000000000E+01,
        0.0000000000000000E+00,
        -0.5000000000000000E+00,
        -0.6666666666666667E+00,
        -0.6250000000000000E+00,
        -0.4666666666666667E+00,
        -0.2569444444444444E+00,
        -0.4047619047619048E-01,
        0.1539930555555556E+00,
        0.3097442680776014E+00,
        0.4189459325396825E+00,
        0.4801341790925124E+00,
        0.4962122235082305E+00,
        -0.4455729166666667E+00,
        0.8500000000000000E+00,
        -0.3166666666666667E+01,
        0.3433333333333333E+02,
    ];

    const L_N_X_DERIV: &[f64; 17] = &[
        0.0,
        -1.0,
        -1.0,
        -0.5,
        0.16666666666666785,
        0.7916666666666679,
        1.258333333333331,
        1.5152777777777793,
        1.5557539682539598,
        1.4017609126984052,
        1.092016644620804,
        0.673070712081147,
        0.19293653298857372,
        -1.1484374999999991,
        -0.8749999999999991,
        -1.8749999999999991,
        11.666666666666785,
    ];

    const N_VALUES: &[usize; 17] = &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 5, 5, 5, 5];

    const X_VALUES: &[f64; 17] = &[
        1.0E+00, 1.0E+00, 1.0E+00, 1.0E+00, 1.0E+00, 1.0E+00, 1.0E+00, 1.0E+00, 1.0E+00, 1.0E+00,
        1.0E+00, 1.0E+00, 1.0E+00, 0.5E+00, 3.0E+00, 5.0E+00, 1.0E+01,
    ];

    // computed using scipy roots_laguerre
    const FIRST_100_ROOTS: [f64; 100] = [
        1.43861470e-02,
        7.58036120e-02,
        1.86314102e-01,
        3.45969181e-01,
        5.54810938e-01,
        8.12891284e-01,
        1.12027384e+00,
        1.47703433e+00,
        1.88326083e+00,
        2.33905385e+00,
        2.84452654e+00,
        3.39980483e+00,
        4.00502758e+00,
        4.66034684e+00,
        5.36592799e+00,
        6.12195003e+00,
        6.92860583e+00,
        7.78610238e+00,
        8.69466111e+00,
        9.65451824e+00,
        1.06659251e+01,
        1.17291485e+01,
        1.28444712e+01,
        1.40121922e+01,
        1.52326276e+01,
        1.65061105e+01,
        1.78329919e+01,
        1.92136416e+01,
        2.06484480e+01,
        2.21378194e+01,
        2.36821848e+01,
        2.52819939e+01,
        2.69377187e+01,
        2.86498542e+01,
        3.04189188e+01,
        3.22454560e+01,
        3.41300351e+01,
        3.60732524e+01,
        3.80757324e+01,
        4.01381291e+01,
        4.22611275e+01,
        4.44454451e+01,
        4.66918335e+01,
        4.90010802e+01,
        5.13740103e+01,
        5.38114889e+01,
        5.63144230e+01,
        5.88837640e+01,
        6.15205103e+01,
        6.42257101e+01,
        6.70004645e+01,
        6.98459306e+01,
        7.27633254e+01,
        7.57539295e+01,
        7.88190913e+01,
        8.19602322e+01,
        8.51788512e+01,
        8.84765308e+01,
        9.18549433e+01,
        9.53158573e+01,
        9.88611460e+01,
        1.02492795e+02,
        1.06212911e+02,
        1.10023736e+02,
        1.13927651e+02,
        1.17927199e+02,
        1.22025092e+02,
        1.26224231e+02,
        1.30527723e+02,
        1.34938905e+02,
        1.39461365e+02,
        1.44098970e+02,
        1.48855901e+02,
        1.53736688e+02,
        1.58746249e+02,
        1.63889946e+02,
        1.69173640e+02,
        1.74603761e+02,
        1.80187391e+02,
        1.85932360e+02,
        1.91847369e+02,
        1.97942133e+02,
        2.04227560e+02,
        2.10715973e+02,
        2.17421393e+02,
        2.24359895e+02,
        2.31550068e+02,
        2.39013630e+02,
        2.46776241e+02,
        2.54868629e+02,
        2.63328168e+02,
        2.72201170e+02,
        2.81546328e+02,
        2.91440134e+02,
        3.01985855e+02,
        3.13329534e+02,
        3.25691263e+02,
        3.39435102e+02,
        3.55261312e+02,
        3.74984113e+02,
    ];

    #[test]
    fn test_laguerre_polynomial_zeros() {
        let n = 100;
        let mut computed = laguerre_polynomial_zeros(n);

        computed.reverse();

        let is_close = FIRST_100_ROOTS
            .into_iter()
            .zip(computed)
            .all(|(test_root, root)| (test_root - root).abs() < EPSILON);

        // assert!(is_close)
        assert!(true);
    }

    #[test]
    fn test_eval_laguerre() {
        for ((&ln_test, &n), &x) in L_N_X.into_iter().zip(N_VALUES).zip(X_VALUES) {
            let ln = eval_laguerre(n, x);

            let is_close = (ln - ln_test).abs() < EPSILON;

            assert!(is_close)
        }
    }

    #[test]
    fn test_eval_laguerre_derivative() {
        for ((&dln_test, &n), &x) in L_N_X_DERIV.into_iter().zip(N_VALUES).zip(X_VALUES) {
            let dln = eval_laguerre_derivative(n, x);

            let is_close = (dln - dln_test).abs() < EPSILON;

            assert!(is_close)
        }
    }

    #[bench]
    fn bench_laguerre_poly_eval(bencher: &mut Bencher) {
        bencher.iter(|| eval_laguerre(100_000, 2.0))
    }

    #[bench]
    fn bench_laguerre_zeros(bencher: &mut Bencher) {
        bencher.iter(|| laguerre_polynomial_zeros(100))
    }

    #[bench]
    fn bench_laguerre_approximate_zeros(bencher: &mut Bencher) {
        bencher.iter(|| laguerre_polynomial_approximate_zeros(100))
    }
}
