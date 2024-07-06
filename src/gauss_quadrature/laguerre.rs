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

use num::{Float, ToPrimitive};

use rayon::iter::{IntoParallelIterator, ParallelIterator};

use super::bessel::bessel_j0_zeros;

fn max<F: Float>(a: F, b: F) -> F {
    let two: F = F::one() + F::one();
    ((a + b) + (a - b).abs()) / two
}

fn min<F: Float>(a: F, b: F) -> F {
    let two: F = F::one() + F::one();
    ((a + b) - (a - b).abs()) / two
}

pub fn nb_eigenvalues_lt_x<F: Float + Send + Sync>(
    diagonal: &[F],
    off_diagonal: &[F],
    x: F,
) -> usize {
    let mut q = F::one();
    let epsilon = F::from(f64::EPSILON).unwrap();
    let mut k: usize = 0;
    let n = diagonal.len();

    for i in 0..n {
        q = if q.is_zero() {
            diagonal[i] - x - off_diagonal[i].abs() / epsilon
        } else {
            diagonal[i] - x - off_diagonal[i].powi(2) / q
        };

        if q.is_sign_negative() {
            k += 1;
        }
    }

    k
}

fn gershgorin_bounds<F: Float + Send + Sync>(diagonal: &[F], off_diagonal: &[F]) -> (F, F) {
    let n = diagonal.len();

    // off_diagonal[0] = zero();

    let (lower_bound, upper_bound) = (0..n - 1)
        .into_par_iter()
        .map(|i| {
            let x = off_diagonal[i].abs() + off_diagonal[i + 1].abs();
            (diagonal[i] - x, diagonal[i] + x)
        })
        .reduce(
            || {
                (
                    diagonal[n - 1] - off_diagonal[n - 1].abs(),
                    diagonal[n - 1] + off_diagonal[n - 1].abs(),
                )
            },
            |(l1, u1), (l2, u2)| (min(l1, l2), max(u1, u2)),
        );

    (lower_bound, upper_bound)
}

fn kth_eigenvalue<F: Float + Send + Sync>(diagonal: &[F], off_diagonal: &[F], k: usize) -> F {
    let n = diagonal.len();
    let epsilon = F::from(f64::EPSILON).unwrap();
    let two = F::one() + F::one();

    let (mut xlower, mut xupper) = gershgorin_bounds(diagonal, off_diagonal);

    let mut tolerance = epsilon * (xupper.abs() + xlower.abs());

    while (xupper - xlower).abs() > tolerance {
        let xmid = (xupper + xlower) / two;

        let nb_eig_lt_xmid = nb_eigenvalues_lt_x(diagonal, off_diagonal, xmid);

        if nb_eig_lt_xmid >= n - k {
            xupper = xmid;
        } else {
            xlower = xmid;
        }

        tolerance = epsilon * (xupper.abs() + xlower.abs());
    }

    (xlower + xupper) / two
}

fn eigenvalues<F: Float + Send + Sync>(diagonal: &[F], off_diagonal: &[F]) -> Vec<F> {
    let n = diagonal.len();
    let eigenvalues: Vec<F> = (0..n)
        .into_par_iter()
        .map(|k| kth_eigenvalue(diagonal, off_diagonal, k))
        .collect();
    eigenvalues
}

pub fn laguerre_polynomial_zeros(n: usize) -> Vec<f64> {
    // define the Jacobi matrix (tridiagonal symmetric matrix)

    // we first define the sub-diagonal
    let off_diagonal: Vec<f64> = (0..n).into_par_iter().map(|i| i as f64).collect();

    // then the diagonal
    let diagonal: Vec<f64> = (0..n)
        .into_par_iter()
        .map(|i| {
            let e = 2 * i + 1;
            e as f64
        })
        .collect();

    let zeros = eigenvalues(diagonal.as_slice(), off_diagonal.as_slice());

    zeros
}

pub fn laguerre_polynomial_approximate_zeros(n: usize) -> Vec<f64> {
    (1..=n)
        .into_par_iter()
        .map(|m| approximate_laguerre_zero(m, n))
        .collect()
}

fn approximate_laguerre_zero(m: usize, n: usize) -> f64 {
    let n_f = n.to_f64().unwrap();

    let j_0_m = bessel_j0_zeros(m);
    let k_n: f64 = n_f + 0.5;

    let term1 = j_0_m.powi(2) / (4.0 * k_n);
    let term2 = 1.0 + (j_0_m.powi(2) - 2.0) / (48.0 * k_n.powi(2));

    term1 * term2
}

fn laguerre_zeros_bounds(m: usize, n: usize) -> (f64, f64) {
    // converting n to f64
    let n_f = n.to_f64().unwrap();

    // converting m to f64
    let m_f = m.to_f64().unwrap();

    // Bessel J0 function m'th positive zero
    let j_0_m = bessel_j0_zeros(m);

    let k_n: f64 = n_f + 0.5;
    let k_m: f64 = m_f + 0.5;

    let lower_bound = j_0_m.powi(2) / (4.0 * k_n);

    // computing second term in formula for simplifocation
    let upper_bound_term = (4.0 * k_m + 0.25).sqrt();
    let upper_bound = (k_m / k_n) * (2.0 * k_m + upper_bound_term);

    (lower_bound, upper_bound)
}

#[cfg(test)]
mod tests {

    use super::*;
    use test::Bencher;

    const EPSILON: f64 = 10e-7;

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

        assert!(is_close)
    }

    #[test]
    fn test_laguerre_polynomial_approximate_zeros() {
        let n = 100;
        // let roots = laguerre_polynomial_zeros(n);
        let roots = laguerre_polynomial_approximate_zeros(n);
        for m in 0..n {
            let root = FIRST_100_ROOTS[m];
            let approximated_root = roots[m];

            println!("approximation: {}, root: {}", approximated_root, root);
        }
    }

    #[bench]
    fn bench_laguerre_zeros(bencher: &mut Bencher) {
        bencher.iter(|| laguerre_polynomial_zeros(10_000))
    }

    #[bench]
    fn bench_laguerre_approximate_zeros(bencher: &mut Bencher) {
        bencher.iter(|| laguerre_polynomial_approximate_zeros(10_000))
    }
}
