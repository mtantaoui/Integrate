use integrate::batch::*;
use integrate::gauss_kronrod::gauss_kronrod_rule;
use integrate::gauss_quadrature::{
    gauss_first_kind_chebyshev_rule, gauss_hermite_rule, gauss_laguerre_rule,
    gauss_second_kind_chebyshev_rule, legendre_rule,
};
use integrate::newton_cotes::{newton_rule, rectangle_rule, simpson_rule, trapezoidal_rule};
use integrate::romberg::romberg_method;
use simdly::SliceMath;

fn exp_batch(xs: &mut [f64]) {
    xs.iter_mut().for_each(|x| *x = x.exp());
}

fn cos_batch(xs: &mut [f64]) {
    xs.iter_mut().for_each(|x| *x = x.cos());
}

/// Batch and scalar rules share one implementation, so for an `f64`
/// integrand computing the same values the results must be bit-identical.
fn assert_close(batch: f64, scalar: f64, ctx: &str) {
    assert_eq!(
        batch.to_bits(),
        scalar.to_bits(),
        "{ctx}: batch={batch}, scalar={scalar}"
    );
}

/// Interval counts whose node counts straddle chunk boundaries for each rule
/// (n+1, n, 2n+1, 3n+1 nodes).
fn sizes() -> Vec<usize> {
    let c = CHUNK_SIZE;
    vec![
        1,
        2,
        7,
        100,
        101,
        c / 3,
        c / 2,
        c - 2,
        c - 1,
        c,
        c + 1,
        3 * c + 5,
    ]
}

#[test]
fn newton_cotes_batch_match_scalar() {
    let (a, b) = (-1.0, 2.0);
    for n in sizes() {
        assert_close(
            rectangle_rule_batch(exp_batch, a, b, n),
            rectangle_rule(|x: f64| x.exp(), a, b, n),
            &format!("rectangle n={n}"),
        );
        assert_close(
            trapezoidal_rule_batch(exp_batch, a, b, n),
            trapezoidal_rule(|x: f64| x.exp(), a, b, n),
            &format!("trapezoidal n={n}"),
        );
        assert_close(
            simpson_rule_batch(exp_batch, a, b, n),
            simpson_rule(|x: f64| x.exp(), a, b, n),
            &format!("simpson n={n}"),
        );
        assert_close(
            newton_rule_batch(exp_batch, a, b, n),
            newton_rule(|x: f64| x.exp(), a, b, n),
            &format!("newton n={n}"),
        );
    }
}

#[test]
fn romberg_batch_matches_scalar() {
    // Level 13 evaluates 4096 new midpoints → exercises chunking.
    for n in [1_usize, 2, 5, 10, 14] {
        assert_close(
            romberg_method_batch(exp_batch, -1.0, 2.0, n),
            romberg_method(|x: f64| x.exp(), -1.0, 2.0, n),
            &format!("romberg n={n}"),
        );
    }
}

#[test]
fn legendre_batch_matches_scalar() {
    for n in sizes() {
        assert_close(
            legendre_rule_batch(exp_batch, -1.0, 2.0, n),
            legendre_rule(|x: f64| x.exp(), -1.0, 2.0, n),
            &format!("legendre n={n}"),
        );
    }
}

#[test]
fn weighted_gauss_batch_match_scalar() {
    for n in [1_usize, 2, 10, 30, 60] {
        assert_close(
            gauss_laguerre_rule_batch(cos_batch, n),
            gauss_laguerre_rule(|x: f64| x.cos(), n),
            &format!("laguerre n={n}"),
        );
        assert_close(
            gauss_hermite_rule_batch(cos_batch, n),
            gauss_hermite_rule(|x: f64| x.cos(), n),
            &format!("hermite n={n}"),
        );
    }
    for n in sizes() {
        assert_close(
            gauss_first_kind_chebyshev_rule_batch(cos_batch, n),
            gauss_first_kind_chebyshev_rule(|x: f64| x.cos(), n),
            &format!("chebyshev1 n={n}"),
        );
        assert_close(
            gauss_second_kind_chebyshev_rule_batch(cos_batch, n),
            gauss_second_kind_chebyshev_rule(|x: f64| x.cos(), n),
            &format!("chebyshev2 n={n}"),
        );
    }
}

#[test]
fn kronrod_batch_matches_scalar() {
    // The scalar Kronrod node computation does not converge for n >= 39.
    for n in [1_usize, 2, 7, 15, 30, 38] {
        let (bi, be) = gauss_kronrod_rule_batch(exp_batch, -1.0, 2.0, n).unwrap();
        let (si, se) = gauss_kronrod_rule(|x: f64| x.exp(), -1.0, 2.0, n).unwrap();
        assert_close(bi, si, &format!("kronrod integral n={n}"));
        assert!((be - se).abs() < 1e-12, "kronrod error n={n}: {be} vs {se}");
    }
}

#[test]
fn batch_calls_respect_chunk_size_and_cover_all_nodes() {
    fn counter(seen: &mut usize) -> impl FnMut(&mut [f64]) + '_ {
        move |xs: &mut [f64]| {
            assert!(!xs.is_empty() && xs.len() <= CHUNK_SIZE);
            *seen += xs.len();
        }
    }

    let n = 3 * CHUNK_SIZE + 5;
    let mut seen = 0;

    rectangle_rule_batch(counter(&mut seen), 0.0, 1.0, n);
    assert_eq!(std::mem::take(&mut seen), n);
    trapezoidal_rule_batch(counter(&mut seen), 0.0, 1.0, n);
    assert_eq!(std::mem::take(&mut seen), n + 1);
    simpson_rule_batch(counter(&mut seen), 0.0, 1.0, n);
    assert_eq!(std::mem::take(&mut seen), 2 * n + 1);
    newton_rule_batch(counter(&mut seen), 0.0, 1.0, n);
    assert_eq!(std::mem::take(&mut seen), 3 * n + 1);
    legendre_rule_batch(counter(&mut seen), 0.0, 1.0, n);
    assert_eq!(std::mem::take(&mut seen), n);
    gauss_first_kind_chebyshev_rule_batch(counter(&mut seen), n);
    assert_eq!(std::mem::take(&mut seen), n);
    gauss_kronrod_rule_batch(counter(&mut seen), 0.0, 1.0, 30).unwrap();
    assert_eq!(std::mem::take(&mut seen), 61);
    // Romberg: 2 endpoints + 1 + 2 + ... + 2^(n-2) midpoints = 2^(n-1) + 1.
    romberg_method_batch(counter(&mut seen), 0.0, 1.0, 14_usize);
    assert_eq!(std::mem::take(&mut seen), (1 << 13) + 1);
}

#[test]
fn simdly_integrands_are_accurate() {
    let exact = std::f64::consts::E - 1.0;
    let simd_exp = |xs: &mut [f64]| xs.exp_in_place();

    let trap = trapezoidal_rule_batch(simd_exp, 0.0, 1.0, 100_000_usize);
    assert!((trap - exact).abs() < 1e-10, "trapezoid: {trap}");

    let simp = simpson_rule_batch(simd_exp, 0.0, 1.0, 1_000_usize);
    assert!((simp - exact).abs() < 1e-13, "simpson: {simp}");

    let newt = newton_rule_batch(simd_exp, 0.0, 1.0, 1_000_usize);
    assert!((newt - exact).abs() < 1e-13, "newton: {newt}");

    let romb = romberg_method_batch(simd_exp, 0.0, 1.0, 10_usize);
    assert!((romb - exact).abs() < 1e-13, "romberg: {romb}");

    let (kr, err) = gauss_kronrod_rule_batch(simd_exp, 0.0, 1.0, 7).unwrap();
    assert!(
        (kr - exact).abs() < 1e-13 && err < 1e-12,
        "kronrod: {kr} ± {err}"
    );

    for n in [10_usize, 20, 100, 500, 10_000] {
        let gl = legendre_rule_batch(simd_exp, 0.0, 1.0, n);
        assert!((gl - exact).abs() < 1e-13, "legendre n={n}: {gl}");
    }

    // ∫₀^π sin x dx = 2
    let pi = std::f64::consts::PI;
    let gl = legendre_rule_batch(|xs: &mut [f64]| xs.sin_in_place(), 0.0, pi, 50_usize);
    assert!((gl - 2.0).abs() < 1e-13, "sin: {gl}");

    // ∫ e^{-x²} cos x dx = √π e^{-1/4}. The Hermite nodes themselves are only
    // accurate to ~1e-9 (scalar rule has the same error), hence the tolerance.

    let herm: f64 = gauss_hermite_rule_batch(|xs: &mut [f64]| xs.cos_in_place(), 15);
    let exact_h = pi.sqrt() * (-0.25_f64).exp();
    assert!((herm - exact_h).abs() < 1e-8, "hermite: {herm}");
}

#[test]
fn f32_batch_integrands() {
    let exact = std::f64::consts::E - 1.0;
    let simd_exp = |xs: &mut [f32]| xs.exp_in_place();

    let gl = legendre_rule_batch(simd_exp, 0.0_f32, 1.0, 10_usize);
    assert!((gl - exact).abs() < 1e-6, "legendre: {gl}");

    let simp = simpson_rule_batch(simd_exp, 0.0_f32, 1.0, 100_usize);
    assert!((simp - exact).abs() < 1e-5, "simpson: {simp}");

    let cheb: f32 = gauss_first_kind_chebyshev_rule_batch(|xs: &mut [f32]| xs.cos_in_place(), 20);
    // ∫ cos x / √(1-x²) dx = π J₀(1)
    assert!((cheb as f64 - std::f64::consts::PI * 0.765_197_686_557_966_6).abs() < 1e-5);
}

#[test]
#[should_panic]
fn trapezoidal_batch_zero_intervals_panics() {
    trapezoidal_rule_batch(exp_batch, 0.0, 1.0, 0_usize);
}

#[test]
#[should_panic]
fn simpson_batch_inverted_limits_panics() {
    simpson_rule_batch(exp_batch, 1.0, 0.0, 10_usize);
}

#[test]
#[should_panic]
fn legendre_batch_inverted_limits_panics() {
    legendre_rule_batch(exp_batch, 1.0, 0.0, 10_usize);
}

#[test]
#[should_panic]
fn laguerre_batch_zero_nodes_panics() {
    let _: f64 = gauss_laguerre_rule_batch(exp_batch, 0);
}

#[test]
fn f32_gauss_rules_accumulate_like_scalar() {
    // Both variants accumulate in f64, so they agree exactly for f32 too.
    for n in [10_usize, 1_000, 100_000] {
        let batch: f32 = gauss_first_kind_chebyshev_rule_batch(
            |xs: &mut [f32]| xs.iter_mut().for_each(|x| *x = x.cos()),
            n,
        );
        let scalar: f32 = gauss_first_kind_chebyshev_rule(|x: f32| x.cos(), n);
        assert_eq!(
            batch.to_bits(),
            scalar.to_bits(),
            "n={n}: {batch} vs {scalar}"
        );
    }
}

#[test]
#[should_panic(expected = "node count overflows usize")]
fn newton_batch_node_count_overflow_panics() {
    newton_rule_batch(exp_batch, 0.0, 1.0, usize::MAX / 2);
}

#[test]
#[should_panic(expected = "node count overflows usize")]
fn simpson_scalar_node_count_overflow_panics() {
    simpson_rule(|x: f64| x, 0.0, 1.0, usize::MAX / 2 + 1);
}

#[test]
#[should_panic(expected = "number of columns must be at most 32")]
fn romberg_batch_too_many_columns_panics() {
    romberg_method_batch(exp_batch, 0.0, 1.0, 70_usize);
}

#[test]
#[should_panic(expected = "number of columns must be at most 32")]
fn romberg_scalar_too_many_columns_panics() {
    romberg_method(|x: f64| x, 0.0, 1.0, 33_usize);
}
