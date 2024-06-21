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
use std::ops::AddAssign;

use num::Float;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

pub fn sturm_sequence<F: Float + Send + Sync + AddAssign>(
    d: Vec<F>,
    off: Vec<F>,
    x: F,
    n: usize,
) -> usize {
    let epsilon = F::from(f64::EPSILON).unwrap();

    let (_, sign_changes_count) = (0..n)
        .into_par_iter()
        .map_init(
            || (F::one(), 0),
            |(mut q, mut k), i| {
                q = if q.is_zero() {
                    d[i] - x - off[i].abs() / epsilon
                } else {
                    d[i] - x - off[i] * off[i] / q
                };
                if q < F::zero() {
                    k += 1;
                }

                (q, k)
            },
        )
        .reduce(|| (F::one(), 0), |(_, k1), (_, k2)| (F::one(), k1 + k2));

    sign_changes_count
}

fn givens_bisection<F: Float>(
    diagonal: Vec<F>,
    off_diagonal: Vec<F>,
    relative_tolerance: F,
    n: usize,
) {
    todo!()
}
