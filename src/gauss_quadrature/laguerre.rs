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

use std::marker::{Send, Sync};
use std::ops::{AddAssign, MulAssign};

use num::{zero, Float};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

pub fn sturm_sequence<F: Float + Send + Sync + AddAssign>(
    d: &[F],
    off: &[F],
    x: F,
    n: usize,
) -> usize {
    let epsilon = F::from(f64::EPSILON).unwrap();

    let mut q = F::one();
    let mut k: usize = 0;

    for i in 0..n {
        if q.is_zero() {
            q = d[i] - x - off[i].abs() / epsilon;
        } else {
            q = d[i] - x - off[i] * off[i] / q;
        }
        if q < zero() {
            k += 1;
        }
    }
    k
}

pub fn givens_bisection<F: Float + Sync + Send + MulAssign + AddAssign>(
    diagonal: &[F],
    off_diagonal: &mut [F],
    mut relative_tolerance: F,
    n: usize,
) -> Vec<F> {
    // Use Gerschgorin's Theorem to Find Upper and Lower Bounds for
    // All Eigenvalues.
    off_diagonal[0] = zero();

    let mut upper_bound = diagonal[n - 1] + off_diagonal[n - 1].abs();
    let mut lower_bound = diagonal[n - 1] - off_diagonal[n - 1].abs();

    for i in (0..=n - 2).rev() {
        let x = off_diagonal[i].abs() + off_diagonal[i + 1].abs();

        upper_bound = if diagonal[i] + x > upper_bound {
            diagonal[i] + x
        } else {
            upper_bound
        };

        lower_bound = if diagonal[i] - x < lower_bound {
            diagonal[i] - x
        } else {
            lower_bound
        }
    }

    // Calculate tolerances
    let mut epsilon = if upper_bound + lower_bound > zero() {
        upper_bound
    } else {
        lower_bound
    };
    epsilon *= F::from(f64::EPSILON).unwrap();

    relative_tolerance = if relative_tolerance < epsilon {
        epsilon
    } else {
        relative_tolerance
    };

    // Initialize Upperbounds and Lowerbounds

    let mut eigenvalues = vec![upper_bound; n];
    let mut lower_bounds = vec![lower_bound; n];

    // Find all eigenvalues from largest to smallest storing
    // from smallest to largest.

    let mut xupper = upper_bound;
    for k in (0..=n - 1).rev() {
        let mut xlower = lower_bound;

        for i in (0..=k).rev() {
            if xlower < lower_bounds[i] {
                xlower = lower_bounds[i];
                break;
            }
        }

        xupper = if xupper > eigenvalues[k] {
            eigenvalues[k]
        } else {
            xupper
        };

        let mut tolerance =
            F::from(2.0).unwrap() * F::from(f64::EPSILON).unwrap() * (xlower.abs() + xupper.abs())
                + relative_tolerance;

        while xupper - xlower > tolerance {
            let xmid = F::from(0.5).unwrap() * (xupper + xlower);

            let j = sturm_sequence(diagonal, off_diagonal, xmid, n) as isize - 1;

            if j < k as isize {
                if j < 0 {
                    lower_bounds[0] = xmid;
                    xlower = xmid;
                } else {
                    xlower = xmid;
                    lower_bounds[j as usize + 1] = xmid;

                    if eigenvalues[j as usize] > xmid {
                        eigenvalues[j as usize] = xmid;
                    }
                }
            } else {
                xupper = xmid;
            }
            tolerance = F::from(2.0).unwrap()
                * F::from(f64::EPSILON).unwrap()
                * (xlower.abs() + xupper.abs())
                + relative_tolerance;
        }
        eigenvalues[k] = F::from(0.5).unwrap() * (xupper + xlower);
    }
    eigenvalues
}

pub fn laguerre_polynomial_zeros<
    F: Float + Send + Sync + MulAssign + AddAssign + std::fmt::Debug,
>(
    n: usize,
) -> Vec<F> {
    // define the jacobi matrix
    let mut bj: Vec<F> = (0..n)
        .into_par_iter()
        .map(|i| F::from(i).unwrap())
        .collect();

    let x: Vec<F> = (0..n)
        .into_par_iter()
        .map(|i| F::from(2 * i + 1).unwrap())
        .collect();

    let relative_tolerance = F::from(1e-10).unwrap();

    let zeros = givens_bisection(x.as_ref(), &mut bj, relative_tolerance, n);

    zeros
}
