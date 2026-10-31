//! Scalar closures vs. batch integrands (plain loop and simdly SIMD).
//!
//! For each rule and size, three variants integrate the same function:
//! * `scalar`       — existing API, `Fn(f64) -> f64` called once per node;
//! * `batch_loop`   — batch API with a plain per-element loop (measures the
//!   batching overhead/benefit on its own);
//! * `batch_simdly` — batch API with `simdly` `*_in_place` SIMD kernels.

use codspeed_criterion_compat::{
    black_box, criterion_group, criterion_main, BenchmarkId, Criterion,
};
use integrate::prelude::*;
use simdly::SliceMath;

fn exp_scalar(x: f64) -> f64 {
    x.exp()
}
fn exp_loop(xs: &mut [f64]) {
    xs.iter_mut().for_each(|x| *x = x.exp());
}
fn exp_simd(xs: &mut [f64]) {
    xs.exp_in_place();
}
fn cos_scalar(x: f64) -> f64 {
    x.cos()
}
fn cos_loop(xs: &mut [f64]) {
    xs.iter_mut().for_each(|x| *x = x.cos());
}
fn cos_simd(xs: &mut [f64]) {
    xs.cos_in_place();
}

type Batch = fn(&mut [f64]);

/// Benchmarks the scalar rule against the batch rule with a plain loop and with simdly.
fn compare<R1, R2>(
    c: &mut Criterion,
    group: &str,
    sizes: &[usize],
    scalar: impl Fn(usize) -> R1,
    batch: impl Fn(Batch, usize) -> R2,
    (loop_fn, simd_fn): (Batch, Batch),
) {
    let mut g = c.benchmark_group(group);
    for &n in sizes {
        g.bench_with_input(BenchmarkId::new("scalar", n), &n, |b, &n| {
            b.iter(|| scalar(black_box(n)))
        });
        g.bench_with_input(BenchmarkId::new("batch_loop", n), &n, |b, &n| {
            b.iter(|| batch(loop_fn, black_box(n)))
        });
        g.bench_with_input(BenchmarkId::new("batch_simdly", n), &n, |b, &n| {
            b.iter(|| batch(simd_fn, black_box(n)))
        });
    }
    g.finish();
}

const EXP: (Batch, Batch) = (exp_loop, exp_simd);
const COS: (Batch, Batch) = (cos_loop, cos_simd);

fn bench_newton_cotes(c: &mut Criterion) {
    let sizes = [1_000, 100_000];
    let (a, b) = (0.0_f64, 1.0_f64);
    compare(
        c,
        "batch/rectangle/exp",
        &sizes,
        |n| rectangle_rule(exp_scalar, a, b, n),
        |f, n| rectangle_rule_batch(f, a, b, n),
        EXP,
    );
    compare(
        c,
        "batch/trapezoidal/exp",
        &sizes,
        |n| trapezoidal_rule(exp_scalar, a, b, n),
        |f, n| trapezoidal_rule_batch(f, a, b, n),
        EXP,
    );
    compare(
        c,
        "batch/simpson/exp",
        &sizes,
        |n| simpson_rule(exp_scalar, a, b, n),
        |f, n| simpson_rule_batch(f, a, b, n),
        EXP,
    );
    compare(
        c,
        "batch/newton_38/exp",
        &sizes,
        |n| newton_rule(exp_scalar, a, b, n),
        |f, n| newton_rule_batch(f, a, b, n),
        EXP,
    );
    compare(
        c,
        "batch/trapezoidal/cos",
        &sizes,
        |n| trapezoidal_rule(cos_scalar, a, b, n),
        |f, n| trapezoidal_rule_batch(f, a, b, n),
        COS,
    );
}

fn bench_romberg(c: &mut Criterion) {
    let (a, b) = (0.0_f64, 1.0_f64);
    compare(
        c,
        "batch/romberg/exp",
        &[10, 18],
        |n| romberg_method(exp_scalar, a, b, n),
        |f, n| romberg_method_batch(f, a, b, n),
        EXP,
    );
}

fn bench_gauss(c: &mut Criterion) {
    let (a, b) = (0.0_f64, 1.0_f64);
    // <= 100 uses tabulated nodes; > 100 computes nodes asymptotically.
    compare(
        c,
        "batch/legendre/exp",
        &[10, 100, 10_000],
        |n| legendre_rule(exp_scalar, a, b, n),
        |f, n| legendre_rule_batch(f, a, b, n),
        EXP,
    );
    compare(
        c,
        "batch/kronrod/exp",
        &[7, 30],
        |n| gauss_kronrod_rule(exp_scalar, a, b, n),
        |f, n| gauss_kronrod_rule_batch(f, a, b, n),
        EXP,
    );
    compare(
        c,
        "batch/chebyshev_first_kind/cos",
        &[100, 10_000],
        |n| gauss_first_kind_chebyshev_rule(cos_scalar, n),
        |f, n| gauss_first_kind_chebyshev_rule_batch::<_, f64>(f, n),
        COS,
    );
    compare(
        c,
        "batch/chebyshev_second_kind/cos",
        &[100, 10_000],
        |n| gauss_second_kind_chebyshev_rule(cos_scalar, n),
        |f, n| gauss_second_kind_chebyshev_rule_batch::<_, f64>(f, n),
        COS,
    );
    compare(
        c,
        "batch/laguerre/cos",
        &[20],
        |n| gauss_laguerre_rule(cos_scalar, n),
        |f, n| gauss_laguerre_rule_batch::<_, f64>(f, n),
        COS,
    );
    compare(
        c,
        "batch/hermite/cos",
        &[20],
        |n| gauss_hermite_rule(cos_scalar, n),
        |f, n| gauss_hermite_rule_batch::<_, f64>(f, n),
        COS,
    );
}

criterion_group!(benches, bench_newton_cotes, bench_romberg, bench_gauss);
criterion_main!(benches);
