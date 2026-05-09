use codspeed_criterion_compat::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use integrate::newton_cotes::{newton_rule, rectangle_rule, simpson_rule, trapezoidal_rule};

fn bench_smooth(c: &mut Criterion) {
    let f = |x: f64| x.exp();
    let a = 0.0_f64;
    let b = 1.0_f64;

    let mut group = c.benchmark_group("newton_cotes/smooth");

    for n in [100_usize, 1_000, 10_000, 100_000] {
        group.bench_with_input(BenchmarkId::new("rectangle", n), &n, |bench, &n| {
            bench.iter(|| rectangle_rule(black_box(f), black_box(a), black_box(b), n));
        });
        group.bench_with_input(BenchmarkId::new("trapezoidal", n), &n, |bench, &n| {
            bench.iter(|| trapezoidal_rule(black_box(f), black_box(a), black_box(b), n));
        });
        group.bench_with_input(BenchmarkId::new("simpson", n), &n, |bench, &n| {
            bench.iter(|| simpson_rule(black_box(f), black_box(a), black_box(b), n));
        });
        group.bench_with_input(BenchmarkId::new("newton_38", n), &n, |bench, &n| {
            bench.iter(|| newton_rule(black_box(f), black_box(a), black_box(b), n));
        });
    }

    group.finish();
}

fn bench_oscillatory(c: &mut Criterion) {
    use std::f64::consts::PI;
    let f = |x: f64| (50.0 * PI * x).sin().powi(2);
    let a = 0.0_f64;
    let b = 1.0_f64;

    let mut group = c.benchmark_group("newton_cotes/oscillatory");

    for n in [1_000_usize, 10_000, 100_000] {
        group.bench_with_input(BenchmarkId::new("rectangle", n), &n, |bench, &n| {
            bench.iter(|| rectangle_rule(black_box(f), black_box(a), black_box(b), n));
        });
        group.bench_with_input(BenchmarkId::new("trapezoidal", n), &n, |bench, &n| {
            bench.iter(|| trapezoidal_rule(black_box(f), black_box(a), black_box(b), n));
        });
        group.bench_with_input(BenchmarkId::new("simpson", n), &n, |bench, &n| {
            bench.iter(|| simpson_rule(black_box(f), black_box(a), black_box(b), n));
        });
        group.bench_with_input(BenchmarkId::new("newton_38", n), &n, |bench, &n| {
            bench.iter(|| newton_rule(black_box(f), black_box(a), black_box(b), n));
        });
    }

    group.finish();
}

criterion_group!(benches, bench_smooth, bench_oscillatory);
criterion_main!(benches);
