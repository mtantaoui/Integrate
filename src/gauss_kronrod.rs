//! Gauss-Kronrod quadrature for numerical integration.
//!
//! The (2N+1)-point Gauss-Kronrod rule extends an N-point Gauss-Legendre rule
//! by adding N+1 optimally placed Kronrod points.  This gives both a
//! high-accuracy integral estimate and a reliable error bound at roughly the
//! cost of a single (2N+1)-point function evaluation pass.
//!
//! The node and weight computation is ported from the Fortran77 original by
//! Robert Piessens and Maria Branders (C translation by John Burkardt, GNU LGPL):
//!
//! > Piessens, Branders, "A Note on the Optimal Addition of Abscissas to
//! > Quadrature Formulas of Gauss and Lobatto", *Mathematics of Computation*,
//! > Vol. 28, No. 125, 1974, pp. 135–139.

use std::fmt;

use num::{Float, ToPrimitive};

use crate::utils::checkers::check_gauss_rule_with_limits;

/// Error returned when the Gauss-Kronrod abscissa Newton iteration fails to
/// converge within 50 steps.
#[derive(Debug, Clone)]
pub struct GaussKronrodError;

impl fmt::Display for GaussKronrodError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Gauss-Kronrod: abscissa Newton iteration did not converge within 50 steps"
        )
    }
}

type KronrodData = (Vec<f64>, Vec<f64>, Vec<f64>);

/// Shared context for abscissa/weight computation — fixed for a given `n`.
struct NodeCtx<'a> {
    n: usize,
    m: usize,
    eps: f64,
    coef2: f64,
    even: bool,
    b: &'a [f64],
}

/// Computes one Kronrod abscissa and its weight via Newton iteration on the
/// Chebyshev expansion of the Stieltjes polynomial.
///
/// Mirrors `abwe1` from Piessens & Branders.
/// Returns `(converged_abscissa, kronrod_weight)`.
fn abwe1(ctx: &NodeCtx, mut x: f64) -> Result<(f64, f64), GaussKronrodError> {
    let &NodeCtx {
        n,
        m,
        eps,
        coef2,
        even,
        b,
    } = ctx;

    // ka == 1 means the next Newton step will be the last.
    let mut ka: i32 = if x == 0.0 { 1 } else { 0 };
    let mut delta;
    let mut fd = 0.0_f64;

    for _ in 0..50 {
        let mut b0 = 0.0_f64;
        let mut b1 = 0.0_f64;
        let mut b2 = b[m];
        let yy = 4.0 * x * x - 2.0;
        let mut d1 = 0.0_f64;

        let (mut ai, mut d2, dif) = if even {
            let ai0 = (2 * m + 1) as f64;
            (ai0, ai0 * b[m], 2.0_f64)
        } else {
            ((m + 1) as f64, 0.0_f64, 1.0_f64)
        };

        // Clenshaw recurrence for the Stieltjes polynomial and its derivative.
        for k in 1..=m {
            ai -= dif;
            let i = m - k + 1; // 1-indexed
            b0 = b1;
            b1 = b2;
            let d0 = d1;
            d1 = d2;
            b2 = yy * b1 - b0 + b[i - 1];
            // For !even the derivative uses b[i] (1-indexed) = b[i] (0-indexed).
            let bi_d = if even { b[i - 1] } else { b[i] };
            d2 = yy * d1 - d0 + ai * bi_d;
        }

        // Evaluate polynomial (f) and its derivative (fd) at x.
        let f = if even {
            fd = d2 + d1;
            x * (b2 - b1)
        } else {
            fd = 4.0 * x * d2;
            0.5 * (b2 - b0)
        };

        // Newton correction.
        delta = f / fd;
        x -= delta;

        if ka == 1 {
            break;
        }
        if delta.abs() <= eps {
            ka = 1;
        }
    }

    if ka != 1 {
        return Err(GaussKronrodError);
    }

    // Weight via Legendre polynomial recurrence at the converged abscissa.
    let mut d0 = 1.0_f64;
    let mut d1 = x;
    let mut ai = 0.0_f64;
    let mut d2 = 1.0_f64; // initialised; overwritten for n >= 2

    for _ in 2..=n {
        ai += 1.0;
        d2 = ((ai + ai + 1.0) * x * d1 - ai * d0) / (ai + 1.0);
        d0 = d1;
        d1 = d2;
    }

    Ok((x, coef2 / (fd * d2)))
}

/// Computes one Gaussian abscissa and the corresponding Gauss and
/// Gauss-Kronrod weights via Newton iteration on the Legendre polynomial.
///
/// Mirrors `abwe2` from Piessens & Branders.
/// Returns `(converged_abscissa, kronrod_weight, gauss_weight)`.
fn abwe2(ctx: &NodeCtx, mut x: f64) -> Result<(f64, f64, f64), GaussKronrodError> {
    let &NodeCtx {
        n,
        m,
        eps,
        coef2,
        even,
        b,
    } = ctx;

    let mut ka: i32 = if x == 0.0 { 1 } else { 0 };
    let mut delta;
    let mut pd2 = 0.0_f64; // P_n'(x) at convergence
    let mut p0_leg = 1.0_f64; // P_{n-1}(x) at convergence

    for _ in 0..50 {
        let mut p0 = 1.0_f64; // P_0(x)
        let mut p1 = x; // P_1(x)
        let mut pd0 = 0.0_f64;
        let mut pd1 = 1.0_f64;
        let mut p2 = 0.0_f64;

        // Special seed to avoid 0/0 in the first Newton step when n == 1.
        if n <= 1 {
            if f64::EPSILON < x.abs() {
                p2 = (3.0 * x * x - 1.0) / 2.0;
                pd2 = 3.0 * x;
            } else {
                p2 = 3.0 * x;
                pd2 = 3.0;
            }
        }

        // Legendre recurrence: build up P_n(x) and P_n'(x).
        let mut ai = 0.0_f64;
        for _ in 2..=n {
            ai += 1.0;
            p2 = ((ai + ai + 1.0) * x * p1 - ai * p0) / (ai + 1.0);
            pd2 = ((ai + ai + 1.0) * (p1 + x * pd1) - ai * pd0) / (ai + 1.0);
            p0 = p1;
            p1 = p2;
            pd0 = pd1;
            pd1 = pd2;
        }

        // After the loop: p0 = P_{n-1}(x), pd2 = P_n'(x).
        p0_leg = p0;

        delta = p2 / pd2;
        x -= delta;

        if ka == 1 {
            break;
        }
        if delta.abs() <= eps {
            ka = 1;
        }
    }

    if ka != 1 {
        return Err(GaussKronrodError);
    }

    // Gauss-Legendre weight: w = 2 / (n * P_n'(x) * P_{n-1}(x)).
    let w_g = 2.0 / (n as f64 * pd2 * p0_leg);

    // Kronrod weight: use the Chebyshev expansion of the Stieltjes polynomial.
    let mut p0 = 0.0_f64;
    let mut p1 = 0.0_f64;
    let mut p2 = b[m];
    let yy = 4.0 * x * x - 2.0;

    for k in 1..=m {
        let i = m - k + 1; // 1-indexed → 0-indexed below
        p0 = p1;
        p1 = p2;
        p2 = yy * p1 - p0 + b[i - 1];
    }

    let w_gk = if even {
        w_g + coef2 / (pd2 * x * (p2 - p1))
    } else {
        w_g + 2.0 * coef2 / (pd2 * (p2 - p0))
    };

    Ok((x, w_gk, w_g))
}

/// Computes the N+1 nonnegative abscissas and weights for the (2N+1)-point
/// Gauss-Kronrod rule on `[-1, +1]`.
///
/// Returns `(x, w_gk, w_g)`:
/// - `x[i]` — abscissas in **decreasing** order; `x[n] = 0` always.
/// - `w_gk[i]` — Gauss-Kronrod weights.
/// - `w_g[i]` — Gauss weights (zero for Kronrod-only points).
fn kronrod_nodes(n: usize, eps: f64) -> Result<KronrodData, GaussKronrodError> {
    let m = n.div_ceil(2);
    let even = 2 * m == n;

    let mut x = vec![0.0_f64; n + 1];
    let mut w_gk = vec![0.0_f64; n + 1];
    let mut w_g = vec![0.0_f64; n + 1];

    // Chebyshev coefficients of the Stieltjes polynomial (b[0..=m])
    // and scratch space for their recurrence (tau[0..m]).
    let mut b = vec![0.0_f64; m + 1];
    let mut tau = vec![0.0_f64; m.max(1)]; // max(1) avoids zero-length alloc

    let an = n as f64;

    // Build Chebyshev coefficients.
    if m > 0 {
        tau[0] = (an + 2.0) / (an + an + 3.0);
        b[m - 1] = tau[0] - 1.0;
        let mut ak = an;

        for l in 1..m {
            ak += 2.0;
            tau[l] = ((ak - 1.0) * ak - an * (an + 1.0)) * (ak + 2.0) * tau[l - 1]
                / (ak * ((ak + 3.0) * (ak + 2.0) - an * (an + 1.0)));
            b[m - l - 1] = tau[l];
            for ll in 1..=l {
                b[m - l - 1] += tau[ll - 1] * b[m - l + ll - 1];
            }
        }
    }
    b[m] = 1.0;

    // coef2 = 2^(2n+1) * (n!)^2 / (2n+1)!
    let mut coef2 = 2.0 / (2 * n + 1) as f64;
    for i in 1..=n {
        coef2 *= 4.0 * i as f64 / (n + i) as f64;
    }

    // Initial approximate abscissas via a trigonometric grid spaced π/(2n+1).
    let theta0 = std::f64::consts::PI / (2.0 * an + 1.0);
    let bb_init = theta0.sin();
    let x1_init = theta0.cos();
    let s = bb_init; // sin(θ₀) — rotate by θ₀ per node
    let c = x1_init; // cos(θ₀)
    let coef = 1.0 - (1.0 - 1.0 / an) / (8.0 * an * an);

    let ctx = NodeCtx {
        n,
        m,
        eps,
        coef2,
        even,
        b: &b,
    };

    let mut xx = coef * x1_init;
    let mut x1 = x1_init;
    let mut bb = bb_init;

    // Main loop: alternate Kronrod (abwe1) and Gauss (abwe2) abscissas.
    let mut k = 1_usize;
    while k <= n {
        // ── Kronrod abscissa at index k-1 ──────────────────────────────────
        let (xx_new, wk) = abwe1(&ctx, xx)?;
        xx = xx_new;
        w_gk[k - 1] = wk;
        w_g[k - 1] = 0.0;
        x[k - 1] = xx;

        // Rotate the trig approximation by π/(2n+1).
        let y = x1;
        x1 = y * c - bb * s;
        bb = y * s + bb * c;

        // Next starting estimate: 0 if this is the last Gauss abscissa (odd n),
        // otherwise advance the trig grid.
        xx = if k == n { 0.0 } else { coef * x1 };

        // ── Gauss abscissa at index k ───────────────────────────────────────
        let (xx_new, wgk, wg) = abwe2(&ctx, xx)?;
        xx = xx_new;
        w_gk[k] = wgk;
        w_g[k] = wg;
        x[k] = xx;

        // Rotate again for the next iteration.
        let y = x1;
        x1 = y * c - bb * s;
        bb = y * s + bb * c;
        xx = coef * x1;

        k += 2;
    }

    // For even n, one final Kronrod-only abscissa sits exactly at the origin.
    if even {
        let (xx_zero, wk) = abwe1(&ctx, 0.0)?;
        w_gk[n] = wk;
        w_g[n] = 0.0;
        x[n] = xx_zero;
    }

    Ok((x, w_gk, w_g))
}

/// Approximates $\int_a^b f(x)\\,dx$ using the $(2N+1)$-point Gauss-Kronrod rule.
///
/// The rule is constructed from an $N$-point Gauss-Legendre rule by adding
/// $N+1$ Kronrod points.  Both estimates are computed in a single pass over
/// $2N+1$ function evaluations, and their difference gives the error bound.
///
/// # Parameters
///
/// * `func` — Integrand $f$; a single-variable function.
/// * `lower_limit` — Lower limit of integration $a$.
/// * `upper_limit` — Upper limit of integration $b$.
/// * `n` — Order of the underlying Gauss-Legendre rule; the Kronrod rule uses
///   $2N+1$ points.  Typical values: 7 (matches scipy's default `quad`),
///   10, 15.  Must be at least 1.
///
/// # Returns
///
/// `Ok((integral, error))` where `error` is `|GK_estimate − G_estimate|`,
/// a reliable indicator of absolute accuracy for smooth integrands.
///
/// # Panics
///
/// Delegates argument validation to [`check_gauss_rule_with_limits`], which panics if
/// `n` is zero, either limit is infinite or NaN, or `lower_limit > upper_limit`.
///
/// # Errors
///
/// Returns `Err(GaussKronrodError)` if the internal Newton iteration for any
/// abscissa does not converge within 50 steps (extremely rare in practice).
///
/// # Examples
///
/// ```
/// use integrate::gauss_kronrod::gauss_kronrod_rule;
///
/// let f = |x: f64| x.exp();
/// let (integral, error) = gauss_kronrod_rule(f, 0.0, 1.0, 7).unwrap();
///
/// let exact = std::f64::consts::E - 1.0;
/// assert!((integral - exact).abs() < 1e-12);
/// assert!(error < 1e-12);
/// ```
///
/// # Resources
///
/// * [QUADPACK](https://en.wikipedia.org/wiki/QUADPACK)
/// * Piessens, Branders, "A Note on the Optimal Addition of Abscissas to
///   Quadrature Formulas of Gauss and Lobatto", *Math. Comp.* 28 (1974).
pub fn gauss_kronrod_rule<Func, F1, F2>(
    func: Func,
    lower_limit: F1,
    upper_limit: F1,
    n: usize,
) -> Result<(f64, f64), GaussKronrodError>
where
    F1: Float + ToPrimitive,
    F2: Float + ToPrimitive,
    Func: Fn(F1) -> F2,
{
    check_gauss_rule_with_limits(lower_limit, upper_limit, n);

    let a = lower_limit.to_f64().unwrap();
    let b = upper_limit.to_f64().unwrap();

    // Compute nodes and weights on [-1, +1].
    let (x, w_gk, w_g) = kronrod_nodes(n, f64::EPSILON)?;

    // Change of variables: ∫_a^b f(x)dx = half_len * ∫_{-1}^{1} f(mid + half_len*t) dt
    let mid = (a + b) / 2.0;
    let half_len = (b - a) / 2.0;

    let f = |t: f64| -> f64 {
        func(F1::from(t).unwrap())
            .to_f64()
            .expect("failed to convert f(x) to f64")
    };

    // x[n] == 0 always; contributes only once (no symmetric counterpart).
    let mut gk_sum = w_gk[n] * f(mid);
    let mut g_sum = w_g[n] * f(mid);

    // All other stored abscissas x[i] > 0 have a symmetric counterpart at -x[i].
    for i in 0..n {
        let xi = x[i]; // positive abscissa in (0, 1]
        let fp = f(mid + half_len * xi);
        let fm = f(mid - half_len * xi);
        gk_sum += w_gk[i] * (fp + fm);
        g_sum += w_g[i] * (fp + fm);
    }

    let integral = half_len * gk_sum;
    let error = (half_len * (gk_sum - g_sum)).abs();

    Ok((integral, error))
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPSILON: f64 = 1e-10;

    #[test]
    fn test_exp() {
        // ∫₀¹ eˣ dx = e - 1
        let (integral, error) = gauss_kronrod_rule(|x: f64| x.exp(), 0.0, 1.0, 7).unwrap();
        let exact = std::f64::consts::E - 1.0;
        assert!(
            (integral - exact).abs() < EPSILON,
            "integral={integral}, exact={exact}"
        );
        assert!(error < EPSILON);
    }

    #[test]
    fn test_polynomial() {
        // ∫₀¹ x² dx = 1/3
        let (integral, _) = gauss_kronrod_rule(|x: f64| x * x, 0.0, 1.0, 7).unwrap();
        assert!((integral - 1.0 / 3.0).abs() < EPSILON);
    }

    #[test]
    fn test_sin() {
        // ∫₀^π sin(x) dx = 2
        let (integral, error) =
            gauss_kronrod_rule(|x: f64| x.sin(), 0.0, std::f64::consts::PI, 7).unwrap();
        assert!((integral - 2.0).abs() < EPSILON, "integral={integral}");
        assert!(error < EPSILON);
    }

    #[test]
    fn test_symmetric_interval() {
        // ∫₋₁¹ x³ dx = 0  (odd integrand)
        let (integral, _) = gauss_kronrod_rule(|x: f64| x.powi(3), -1.0, 1.0, 7).unwrap();
        assert!(integral.abs() < EPSILON);
    }

    #[test]
    fn test_f32_input() {
        // Verify the generic float path works for f32 input / f64 output.
        let (integral, _) =
            gauss_kronrod_rule(|x: f32| (x * x) as f64, 0.0_f32, 1.0_f32, 7).unwrap();
        assert!((integral - 1.0 / 3.0).abs() < 1e-6);
    }

    #[test]
    fn test_even_n() {
        // n=4 (even) exercises the extra abwe1 call for the origin.
        let (integral, _) = gauss_kronrod_rule(|x: f64| x.exp(), 0.0, 1.0, 4).unwrap();
        let exact = std::f64::consts::E - 1.0;
        assert!((integral - exact).abs() < EPSILON);
    }

    #[test]
    #[should_panic(expected = "number of steps can't be zero")]
    fn test_n_zero() {
        let _ = gauss_kronrod_rule(|x: f64| x, 0.0, 1.0, 0);
    }

    #[test]
    #[should_panic(expected = "integration limits a and b can't be infinite")]
    fn test_infinite_limit() {
        let _ = gauss_kronrod_rule(|x: f64| x, f64::NEG_INFINITY, 1.0, 7);
    }

    #[test]
    #[should_panic(expected = "lower_limit must not exceed upper_limit")]
    fn test_inverted_limits() {
        let _ = gauss_kronrod_rule(|x: f64| x, 1.0, 0.0, 7);
    }
}
