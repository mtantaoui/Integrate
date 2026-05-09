use codspeed_criterion_compat::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use integrate::adaptive_quadrature::adaptive_simpson_method;
use integrate::romberg::romberg_method;

// Compare adaptive Simpson and Romberg on a smooth integrand, sweeping
// over accuracy targets (adaptive) and extrapolation depth (Romberg).
fn bench_smooth(c: &mut Criterion) {
    let f = |x: f64| x.exp();
    let a = 0.0_f64;
    let b = 1.0_f64;
    let min_h = 1e-10_f64;

    let mut group = c.benchmark_group("adaptive/smooth");

    for tol_exp in [4_i32, 6, 8, 10] {
        let tol = 10.0_f64.powi(-tol_exp);
        group.bench_with_input(
            BenchmarkId::new("adaptive_simpson", format!("1e-{tol_exp}")),
            &tol,
            |bench, &tol| {
                bench.iter(|| {
                    adaptive_simpson_method(
                        black_box(f),
                        black_box(a),
                        black_box(b),
                        min_h,
                        tol,
                    )
                    .unwrap()
                });
            },
        );
    }

    // n_columns controls Richardson extrapolation depth; total evaluations ≈ 2^(n-1)+1.
    for n_cols in [4_usize, 6, 8, 10, 12] {
        group.bench_with_input(
            BenchmarkId::new("romberg", n_cols),
            &n_cols,
            |bench, &n| {
                bench.iter(|| {
                    romberg_method::<_, f64, f64, usize>(
                        black_box(f),
                        black_box(a),
                        black_box(b),
                        n,
                    )
                });
            },
        );
    }

    group.finish();
}

// Same comparison on a mildly oscillatory integrand.
fn bench_oscillatory(c: &mut Criterion) {
    use std::f64::consts::PI;
    let f = |x: f64| (10.0 * PI * x).sin();
    let a = 0.0_f64;
    let b = 1.0_f64;
    let min_h = 1e-10_f64;

    let mut group = c.benchmark_group("adaptive/oscillatory");

    for tol_exp in [4_i32, 6, 8] {
        let tol = 10.0_f64.powi(-tol_exp);
        group.bench_with_input(
            BenchmarkId::new("adaptive_simpson", format!("1e-{tol_exp}")),
            &tol,
            |bench, &tol| {
                bench.iter(|| {
                    adaptive_simpson_method(
                        black_box(f),
                        black_box(a),
                        black_box(b),
                        min_h,
                        tol,
                    )
                    .unwrap()
                });
            },
        );
    }

    for n_cols in [4_usize, 6, 8, 10, 12] {
        group.bench_with_input(
            BenchmarkId::new("romberg", n_cols),
            &n_cols,
            |bench, &n| {
                bench.iter(|| {
                    romberg_method::<_, f64, f64, usize>(
                        black_box(f),
                        black_box(a),
                        black_box(b),
                        n,
                    )
                });
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_smooth, bench_oscillatory);
criterion_main!(benches);
