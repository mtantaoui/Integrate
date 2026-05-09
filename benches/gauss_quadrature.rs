use codspeed_criterion_compat::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use integrate::gauss_kronrod::gauss_kronrod_rule;
use integrate::gauss_quadrature::{
    gauss_first_kind_chebyshev_rule, gauss_second_kind_chebyshev_rule, legendre_rule,
};

// Compare Gauss-Legendre and Gauss-Kronrod on a smooth finite-interval integrand.
fn bench_finite_interval(c: &mut Criterion) {
    let f = |x: f64| x.exp();
    let a = 0.0_f64;
    let b = 1.0_f64;

    let mut group = c.benchmark_group("gauss/finite_interval");

    for n in [10_usize, 50, 100, 500] {
        group.bench_with_input(BenchmarkId::new("legendre", n), &n, |bench, &n| {
            bench.iter(|| legendre_rule(black_box(f), black_box(a), black_box(b), n));
        });
    }

    // Gauss-Kronrod: n is the GL base-rule order; 2n+1 total points.
    for n in [7_usize, 10, 15, 20] {
        group.bench_with_input(BenchmarkId::new("kronrod", n), &n, |bench, &n| {
            bench.iter(|| {
                gauss_kronrod_rule(black_box(f), black_box(a), black_box(b), n).unwrap()
            });
        });
    }

    group.finish();
}

// Compare first- and second-kind Chebyshev rules.
// The integrand is the unweighted part f(x); the weight is built into the rule.
fn bench_chebyshev(c: &mut Criterion) {
    let f = |x: f64| x.cos();

    let mut group = c.benchmark_group("gauss/chebyshev");

    for n in [10_usize, 50, 100, 500] {
        group.bench_with_input(BenchmarkId::new("first_kind", n), &n, |bench, &n| {
            bench.iter(|| gauss_first_kind_chebyshev_rule(black_box(f), n));
        });
        group.bench_with_input(BenchmarkId::new("second_kind", n), &n, |bench, &n| {
            bench.iter(|| gauss_second_kind_chebyshev_rule(black_box(f), n));
        });
    }

    group.finish();
}

criterion_group!(benches, bench_finite_interval, bench_chebyshev);
criterion_main!(benches);
