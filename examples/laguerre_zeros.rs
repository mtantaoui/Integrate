use num::Float;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

fn max<F: Float>(a: F, b: F) -> F {
    let two: F = F::one() + F::one();
    ((a + b) + (a - b).abs()) / two
}

fn min<F: Float>(a: F, b: F) -> F {
    let two: F = F::one() + F::one();
    ((a + b) - (a - b).abs()) / two
}

#[time_graph::instrument]
pub fn nb_eigenvalues_lt_x<F: Float + Send + Sync>(
    diagonal: &[F],
    off_diagonal: &[F],
    x: F,
) -> usize {
    let mut q = F::one();
    let epsilon = F::from(f64::EPSILON).unwrap();
    let mut k: usize = 0;
    let n = diagonal.len();

    for i in 0..n {
        q = if q.is_zero() {
            diagonal[i] - x - off_diagonal[i].abs() / epsilon
        } else {
            diagonal[i] - x - off_diagonal[i].powi(2) / q
        };

        if q.is_sign_negative() {
            k += 1;
        }
    }

    k
}

#[time_graph::instrument]
fn gershgorin_bounds<F: Float + Send + Sync>(diagonal: &[F], off_diagonal: &[F]) -> (F, F) {
    let n = diagonal.len();

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
            },
            |(l1, u1), (l2, u2)| (min(l1, l2), max(u1, u2)),
        );

    (lower_bound, upper_bound)
}

#[time_graph::instrument]
fn kth_eigenvalue<F: Float + Send + Sync>(diagonal: &[F], off_diagonal: &[F], k: usize) -> F {
    let n = diagonal.len();
    let epsilon = F::from(f64::EPSILON).unwrap();
    let two = F::one() + F::one();

    let (mut xlower, mut xupper) = gershgorin_bounds(diagonal, off_diagonal);

    let mut tolerance = two * epsilon * (xupper.abs() + xlower.abs());

    while (xupper - xlower).abs() > tolerance {
        let xmid = (xupper + xlower) / two;

        let nb_eig_lt_xmid = nb_eigenvalues_lt_x(diagonal, off_diagonal, xmid);

        if nb_eig_lt_xmid >= n - k {
            xupper = xmid;
        } else {
            xlower = xmid;
        }

        tolerance = epsilon * (xupper.abs() + xlower.abs());
    }

    (xlower + xupper) / two
}

#[time_graph::instrument]
fn eigenvalues<F: Float + Send + Sync>(diagonal: &[F], off_diagonal: &[F]) -> Vec<F> {
    let n = diagonal.len();
    let eigenvalues: Vec<F> = (0..n)
        .into_par_iter()
        .map(|k| kth_eigenvalue(diagonal, off_diagonal, k))
        .collect();
    eigenvalues
}

#[time_graph::instrument]
pub fn laguerre_polynomial_zeros(n: usize) -> Vec<f64> {
    // define the Jacobi matrix (tridiagonal symmetric matrix)

    // we first define the sub-diagonal
    let off_diagonal: Vec<f64> = (0..n).into_par_iter().map(|i| i as f64).collect();

    // then the diagonal
    let diagonal: Vec<f64> = (0..n)
        .into_par_iter()
        .map(|i| {
            let e = 2 * i + 1;
            e as f64
        })
        .collect();

    let zeros = eigenvalues(diagonal.as_slice(), off_diagonal.as_slice());

    return zeros;
}

#[time_graph::instrument]
fn zeros() {
    let n = 1_000;
    laguerre_polynomial_zeros(n);
}

fn main() {
    time_graph::enable_data_collection(true);

    zeros();

    let graph = time_graph::get_full_graph();

    // println!("{}", graph.as_dot());

    // println!("{}", graph.as_json());

    // println!("{}", graph.as_table());

    println!("{}", graph.as_short_table());
}
