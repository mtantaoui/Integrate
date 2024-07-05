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
//!

use num::{zero, Float, One, Zero};
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};
use rayon::slice::ParallelSlice;

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
