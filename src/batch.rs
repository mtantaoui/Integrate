//! Batch-integrand variants of the quadrature rules.
//!
//! The rules in [`newton_cotes`](crate::newton_cotes),
//! [`gauss_quadrature`](crate::gauss_quadrature),
//! [`gauss_kronrod`](crate::gauss_kronrod) and [`romberg`](crate::romberg)
//! call the integrand once per node. The functions here instead hand the
//! integrand a **slice of nodes** and expect it to overwrite each node $x$
//! with $f(x)$ in place. That lets the integrand evaluate many points at once,
//! e.g. with a SIMD math library such as [`simdly`](https://docs.rs/simdly):
//!
//! ```
//! use integrate::batch::trapezoidal_rule_batch;
//! use simdly::SliceMath;
//!
//! // ∫₀¹ eˣ dx = e - 1
//! let integral = trapezoidal_rule_batch(|xs: &mut [f64]| xs.exp_in_place(), 0.0, 1.0, 10_000_usize);
//! assert!((integral - (std::f64::consts::E - 1.0)).abs() < 1e-8);
//! ```
//!
//! Each `_batch` function shares its implementation with its scalar
//! counterpart: same nodes, weights, summation order and `f64` accumulation.
//! For `f64` integrands that compute the same values, the results are
//! bit-for-bit identical. The adaptive Simpson method has
//! no batch variant: it evaluates two new points per refinement step, so there
//! is nothing to batch.
//!
//! # Integrand contract
//!
//! The integrand has the signature `FnMut(&mut [F])`. On entry each element
//! is a node $x_i$; on return it must hold $f(x_i)$. The slice must not be
//! resized or reordered. The integrand may be called **several times** per
//! integral, with slices of arbitrary (non-zero) length: nodes are processed
//! in chunks of at most [`CHUNK_SIZE`] so that memory use stays bounded and
//! each chunk stays in cache.

use std::{fmt::Debug, ops::AddAssign};

use num::{bigint::ToBigInt, Float, ToPrimitive, Unsigned};

use crate::{
    gauss_kronrod::{gauss_kronrod, GaussKronrodError},
    gauss_quadrature::{gauss_rule, legendre},
    newton_cotes::{newton, rectangle, simpson, trapezoidal},
    romberg::romberg,
    utils::{
        chebyshev::{roots_first_kind_chebyshev, roots_second_kind_chebyshev},
        hermite::roots_hermite,
        integrand::Batch,
        laguerre::roots_laguerre,
    },
};

pub use crate::utils::integrand::CHUNK_SIZE;

// ---------------------------------------------------------------------------
// Newton-Cotes
// ---------------------------------------------------------------------------

/// Batch-integrand version of
/// [`rectangle_rule`](crate::newton_cotes::rectangle_rule) (composite midpoint rule).
///
/// See the [module docs](self) for the integrand contract.
///
/// # Panics
/// Panics if `n_intervals` is zero, either limit is infinite or NaN, or
/// `lower_limit > upper_limit`.
///
/// # Examples
/// ```
/// use integrate::batch::rectangle_rule_batch;
///
/// let square = |xs: &mut [f64]| xs.iter_mut().for_each(|x| *x = *x * *x);
/// let integral = rectangle_rule_batch(square, 0.0, 1.0, 1_000_usize);
/// assert!((integral - 1.0 / 3.0).abs() < 1e-6);
/// ```
pub fn rectangle_rule_batch<Func, F, U>(
    func: Func,
    lower_limit: F,
    upper_limit: F,
    n_intervals: U,
) -> f64
where
    Func: FnMut(&mut [F]),
    F: Float,
    U: Unsigned + ToPrimitive + Copy,
{
    rectangle(&mut Batch::new(func), lower_limit, upper_limit, n_intervals)
}

/// Batch-integrand version of
/// [`trapezoidal_rule`](crate::newton_cotes::trapezoidal_rule).
///
/// Approximates $\int_a^b f(x)\\,dx$ with the composite trapezoidal rule
///
/// $$T_h(f) = h\left(\tfrac{1}{2}f(a) + \sum_{i=1}^{n-1} f(a + ih) + \tfrac{1}{2}f(b)\right), \quad h = \frac{b-a}{n}.$$
///
/// See the [module docs](self) for the integrand contract.
///
/// # Panics
/// Panics if `n_intervals` is zero, either limit is infinite or NaN, or
/// `lower_limit > upper_limit`.
///
/// # Examples
/// ```
/// use integrate::batch::trapezoidal_rule_batch;
///
/// let square = |xs: &mut [f64]| xs.iter_mut().for_each(|x| *x = *x * *x);
/// let integral = trapezoidal_rule_batch(square, 0.0, 1.0, 1_000_usize);
/// assert!((integral - 1.0 / 3.0).abs() < 1e-6);
/// ```
pub fn trapezoidal_rule_batch<Func, F, U>(
    func: Func,
    lower_limit: F,
    upper_limit: F,
    n_intervals: U,
) -> f64
where
    Func: FnMut(&mut [F]),
    F: Float,
    U: Unsigned + ToPrimitive + Copy,
{
    trapezoidal(&mut Batch::new(func), lower_limit, upper_limit, n_intervals)
}

/// Batch-integrand version of
/// [`simpson_rule`](crate::newton_cotes::simpson_rule) (composite Simpson's 1/3 rule).
///
/// Uses the nodes $a + jh/2$, $j = 0, \dots, 2n$, with weights
/// $1, 4, 2, 4, \dots, 2, 4, 1$ scaled by $h/6$.
///
/// See the [module docs](self) for the integrand contract.
///
/// # Panics
/// Panics if `n` is zero, either limit is infinite or NaN, or `a > b`.
///
/// # Examples
/// ```
/// use integrate::batch::simpson_rule_batch;
///
/// let cube = |xs: &mut [f64]| xs.iter_mut().for_each(|x| *x = *x * *x * *x);
/// let integral = simpson_rule_batch(cube, 0.0, 1.0, 10_usize);
/// assert!((integral - 0.25).abs() < 1e-12);
/// ```
pub fn simpson_rule_batch<Func, F, U>(func: Func, a: F, b: F, n: U) -> f64
where
    Func: FnMut(&mut [F]),
    F: Float,
    U: Unsigned + ToPrimitive + Copy,
{
    simpson(&mut Batch::new(func), a, b, n)
}

/// Batch-integrand version of
/// [`newton_rule`](crate::newton_cotes::newton_rule) (composite Newton's 3/8 rule).
///
/// Uses the nodes $a + jh/3$, $j = 0, \dots, 3n$, with weights
/// $1, 3, 3, 2, 3, 3, 2, \dots, 3, 3, 1$ scaled by $h/8$.
///
/// See the [module docs](self) for the integrand contract.
///
/// # Panics
/// Panics if `n_intervals` is zero, either limit is infinite or NaN, or
/// `lower_limit > upper_limit`.
///
/// # Examples
/// ```
/// use integrate::batch::newton_rule_batch;
///
/// let cube = |xs: &mut [f64]| xs.iter_mut().for_each(|x| *x = *x * *x * *x);
/// let integral = newton_rule_batch(cube, 0.0, 1.0, 10_usize);
/// assert!((integral - 0.25).abs() < 1e-12);
/// ```
pub fn newton_rule_batch<Func, F, U>(
    func: Func,
    lower_limit: F,
    upper_limit: F,
    n_intervals: U,
) -> f64
where
    Func: FnMut(&mut [F]),
    F: Float,
    U: Unsigned + ToPrimitive + Copy,
{
    newton(&mut Batch::new(func), lower_limit, upper_limit, n_intervals)
}

// ---------------------------------------------------------------------------
// Romberg
// ---------------------------------------------------------------------------

/// Batch-integrand version of [`romberg_method`](crate::romberg::romberg_method).
///
/// Each refinement level evaluates its $2^{i-1}$ new midpoints as one batch
/// (split into chunks of at most [`CHUNK_SIZE`]).
///
/// See the [module docs](self) for the integrand contract.
///
/// # Panics
/// Panics if `n_columns` is zero or exceeds
/// [`MAX_ROMBERG_COLUMNS`](crate::romberg::MAX_ROMBERG_COLUMNS), either limit
/// is infinite or NaN, or `lower_limit > upper_limit`.
///
/// # Examples
/// ```
/// use integrate::batch::romberg_method_batch;
///
/// let square = |xs: &mut [f64]| xs.iter_mut().for_each(|x| *x = *x * *x);
/// let integral = romberg_method_batch(square, 0.0, 1.0, 10_usize);
/// assert!((integral - 1.0 / 3.0).abs() < 1e-12);
/// ```
pub fn romberg_method_batch<Func, F, U>(
    func: Func,
    lower_limit: F,
    upper_limit: F,
    n_columns: U,
) -> f64
where
    Func: FnMut(&mut [F]),
    F: Float,
    U: Unsigned + ToPrimitive + Copy,
{
    romberg(&mut Batch::new(func), lower_limit, upper_limit, n_columns)
}

// ---------------------------------------------------------------------------
// Gaussian quadrature
// ---------------------------------------------------------------------------

/// Batch-integrand version of
/// [`legendre_rule`](crate::gauss_quadrature::legendre_rule).
///
/// Approximates $\int_a^b f(x)\\,dx$ with the $n$-point Gauss-Legendre rule
///
/// $$\int_a^b f(x)\\,dx \approx c\sum_{k=1}^{n} w_k\\, f(c\\,x_k + d), \quad c = \frac{b-a}{2},\\; d = \frac{b+a}{2}.$$
///
/// See the [module docs](self) for the integrand contract.
///
/// # Panics
/// Panics if `n` is zero, either limit is infinite or NaN, or
/// `lower_limit > upper_limit`.
///
/// # Examples
/// ```
/// use integrate::batch::legendre_rule_batch;
///
/// let square = |xs: &mut [f64]| xs.iter_mut().for_each(|x| *x = *x * *x);
/// let integral = legendre_rule_batch(square, 0.0, 1.0, 10_usize);
/// assert!((integral - 1.0 / 3.0).abs() < 1e-12);
/// ```
pub fn legendre_rule_batch<Func, F, U>(func: Func, lower_limit: F, upper_limit: F, n: U) -> f64
where
    Func: FnMut(&mut [F]),
    F: Float,
    U: Unsigned + ToPrimitive + Copy,
{
    legendre(&mut Batch::new(func), lower_limit, upper_limit, n)
}

/// Batch-integrand version of
/// [`gauss_laguerre_rule`](crate::gauss_quadrature::gauss_laguerre_rule):
/// approximates $\int_0^{\infty} e^{-x} f(x)\\,dx$.
///
/// See the [module docs](self) for the integrand contract.
///
/// # Panics
/// Panics if `n` is zero.
///
/// # Examples
/// ```
/// use integrate::batch::gauss_laguerre_rule_batch;
///
/// // ∫₀^∞ e^{-x} x² dx = 2
/// let square = |xs: &mut [f64]| xs.iter_mut().for_each(|x| *x = *x * *x);
/// let integral = gauss_laguerre_rule_batch(square, 10);
/// assert!((integral - 2.0).abs() < 1e-6);
/// ```
pub fn gauss_laguerre_rule_batch<Func, F>(func: Func, n: usize) -> F
where
    Func: FnMut(&mut [F]),
    F: Float + Debug + Sync + Send + AddAssign,
{
    gauss_rule(&mut Batch::new(func), n, roots_laguerre::<F>)
}

/// Batch-integrand version of
/// [`gauss_hermite_rule`](crate::gauss_quadrature::gauss_hermite_rule):
/// approximates $\int_{-\infty}^{\infty} e^{-x^2} f(x)\\,dx$.
///
/// See the [module docs](self) for the integrand contract.
///
/// # Panics
/// Panics if `n` is zero.
///
/// # Examples
/// ```
/// use integrate::batch::gauss_hermite_rule_batch;
///
/// // ∫ e^{-x²} x² dx = √π / 2
/// let square = |xs: &mut [f64]| xs.iter_mut().for_each(|x| *x = *x * *x);
/// let integral = gauss_hermite_rule_batch(square, 10);
/// assert!((integral - std::f64::consts::PI.sqrt() / 2.0).abs() < 1e-6);
/// ```
pub fn gauss_hermite_rule_batch<Func, F>(func: Func, n: usize) -> F
where
    Func: FnMut(&mut [F]),
    F: Float + Debug + Sync + Send + AddAssign + ToBigInt,
{
    gauss_rule(&mut Batch::new(func), n, roots_hermite::<F>)
}

/// Batch-integrand version of
/// [`gauss_first_kind_chebyshev_rule`](crate::gauss_quadrature::gauss_first_kind_chebyshev_rule):
/// approximates $\int_{-1}^{1} f(x) / \sqrt{1-x^2}\\,dx$.
///
/// See the [module docs](self) for the integrand contract.
///
/// # Panics
/// Panics if `n` is zero.
///
/// # Examples
/// ```
/// use integrate::batch::gauss_first_kind_chebyshev_rule_batch;
///
/// // ∫ x² / √(1-x²) dx = π / 2
/// let square = |xs: &mut [f64]| xs.iter_mut().for_each(|x| *x = *x * *x);
/// let integral = gauss_first_kind_chebyshev_rule_batch(square, 10);
/// assert!((integral - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
/// ```
pub fn gauss_first_kind_chebyshev_rule_batch<Func, F>(func: Func, n: usize) -> F
where
    Func: FnMut(&mut [F]),
    F: Float + Debug + Sync + Send + AddAssign,
{
    gauss_rule(&mut Batch::new(func), n, roots_first_kind_chebyshev::<F>)
}

/// Batch-integrand version of
/// [`gauss_second_kind_chebyshev_rule`](crate::gauss_quadrature::gauss_second_kind_chebyshev_rule):
/// approximates $\int_{-1}^{1} f(x) \sqrt{1-x^2}\\,dx$.
///
/// See the [module docs](self) for the integrand contract.
///
/// # Panics
/// Panics if `n` is zero.
///
/// # Examples
/// ```
/// use integrate::batch::gauss_second_kind_chebyshev_rule_batch;
///
/// // ∫ x² √(1-x²) dx = π / 8
/// let square = |xs: &mut [f64]| xs.iter_mut().for_each(|x| *x = *x * *x);
/// let integral = gauss_second_kind_chebyshev_rule_batch(square, 10);
/// assert!((integral - std::f64::consts::PI / 8.0).abs() < 1e-12);
/// ```
pub fn gauss_second_kind_chebyshev_rule_batch<Func, F>(func: Func, n: usize) -> F
where
    Func: FnMut(&mut [F]),
    F: Float + Debug + Sync + Send + AddAssign,
{
    gauss_rule(&mut Batch::new(func), n, roots_second_kind_chebyshev::<F>)
}

// ---------------------------------------------------------------------------
// Gauss-Kronrod
// ---------------------------------------------------------------------------

/// Batch-integrand version of
/// [`gauss_kronrod_rule`](crate::gauss_kronrod::gauss_kronrod_rule).
///
/// Evaluates the integrand on all $2n+1$ Kronrod nodes and returns
/// `(integral, error_estimate)`, where the error estimate is the difference
/// between the Kronrod and embedded Gauss results.
///
/// See the [module docs](self) for the integrand contract.
///
/// # Errors
/// Returns [`GaussKronrodError`] if the abscissa computation does not converge.
///
/// # Panics
/// Panics if `n` is zero, either limit is infinite or NaN, or
/// `lower_limit > upper_limit`.
///
/// # Examples
/// ```
/// use integrate::batch::gauss_kronrod_rule_batch;
///
/// let exp = |xs: &mut [f64]| xs.iter_mut().for_each(|x| *x = x.exp());
/// let (integral, error) = gauss_kronrod_rule_batch(exp, 0.0, 1.0, 7).unwrap();
/// assert!((integral - (std::f64::consts::E - 1.0)).abs() < 1e-12);
/// assert!(error < 1e-12);
/// ```
pub fn gauss_kronrod_rule_batch<Func, F>(
    func: Func,
    lower_limit: F,
    upper_limit: F,
    n: usize,
) -> Result<(f64, f64), GaussKronrodError>
where
    Func: FnMut(&mut [F]),
    F: Float,
{
    gauss_kronrod(&mut Batch::new(func), lower_limit, upper_limit, n)
}
