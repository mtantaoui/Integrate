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

/// Computes Sturm element associatd with characteristic polynomial
/// using recursive formula.
fn sturm_element<F: Float + Send + Sync + std::fmt::Debug>(
    diagonal: &[F],
    off_diagonal: &[F],
    x: F,
    n: usize,
) -> F {
    if n.is_zero() {
        return F::one();
    }

    if n.is_one() {
        return diagonal[0] - x;
    }
    let (p_r_1, p_r_2) = rayon::join(
        || sturm_element(diagonal, off_diagonal, x, n - 1),
        || sturm_element(diagonal, off_diagonal, x, n - 2),
    );

    (diagonal[n - 1] - x) * p_r_1 - off_diagonal[n - 1].powi(2) * p_r_2
}

pub fn sturm_sequence<F: Float + Send + Sync + std::fmt::Debug>(
    diagonal: &[F],
    off_diagonal: &[F],
    x: F,
) -> Vec<F> {
    let mut sequence = Vec::new();
    let n = diagonal.len();

    (0..n + 1)
        .into_par_iter()
        .map(|i| sturm_element(diagonal, off_diagonal, x, i))
        .collect_into_vec(&mut sequence);

    sequence
}

pub fn nb_eigenvalues_lt_x<F: Float + Send + Sync + std::fmt::Debug>(
    diagonal: &[F],
    off_diagonal: &mut [F],
    x: F,
) -> usize {
    let sturm_seq = sturm_sequence(diagonal, off_diagonal, x);

    let nb_sign_changes: usize = sturm_seq
        .par_windows(2)
        .map(|window| if window[0] * window[1] < zero() { 1 } else { 0 })
        .sum();

    let (lower_bound, upper_bound) = gershgorin_bounds(diagonal, off_diagonal);

    println!("lower : {:?} \nupper : {:?}", lower_bound, upper_bound);

    nb_sign_changes
}

fn gershgorin_bounds<F: Float + Send + Sync>(diagonal: &[F], off_diagonal: &mut [F]) -> (F, F) {
    let n = diagonal.len();

    off_diagonal[0] = zero();

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
            }, // the "identity" is 0 in both columns
            |(l1, u1), (l2, u2)| (min(l1, l2), max(u1, u2)),
        );

    (lower_bound, upper_bound)
}
