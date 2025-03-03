use num::Float;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

// Computes the maximum of two floating point values without using conditionals
// Uses the mathematical identity max(a,b) = (a+b+|a-b|)/2
fn max<F: Float>(a: F, b: F) -> F {
    let two: F = F::one() + F::one();
    ((a + b) + (a - b).abs()) / two
}

// Computes the minimum of two floating point values without using conditionals
// Uses the mathematical identity min(a,b) = (a+b-|a-b|)/2
fn min<F: Float>(a: F, b: F) -> F {
    let two: F = F::one() + F::one();
    ((a + b) - (a - b).abs()) / two
}

// Counts the number of eigenvalues less than a given value x
// Uses a numerical technique based on Sturm sequences
// This is an implementation of the Sturm count for a tridiagonal matrix
#[time_graph::instrument]
pub fn nb_eigenvalues_lt_x<F: Float + Send + Sync>(
    diagonal: &[F],     // The diagonal elements of the tridiagonal matrix
    off_diagonal: &[F], // The off-diagonal elements
    x: F,               // The value to compare eigenvalues against
) -> usize {
    let mut q = F::one();
    let epsilon = F::from(f64::EPSILON).unwrap();
    let mut k: usize = 0; // Counter for eigenvalues less than x
    let n = diagonal.len();

    for i in 0..n {
        // Calculate the next value in the Sturm sequence
        // This handles the special case when q is zero to avoid division by zero
        q = if q.is_zero() {
            diagonal[i] - x - off_diagonal[i].abs() / epsilon
        } else {
            diagonal[i] - x - off_diagonal[i].powi(2) / q
        };

        // Each sign change in the Sturm sequence corresponds to an eigenvalue less than x
        if q.is_sign_negative() {
            k += 1;
        }
    }
    k
}

// Calculates the Gershgorin bounds for the eigenvalues of a tridiagonal matrix
// These bounds provide guaranteed intervals containing all eigenvalues
#[time_graph::instrument]
fn gershgorin_bounds<F: Float + Send + Sync>(diagonal: &[F], off_diagonal: &[F]) -> (F, F) {
    let n = diagonal.len();

    // Compute the bounds in parallel for all but the last element
    let (lower_bound, upper_bound) = (0..n - 1)
        .into_par_iter()
        .map(|i| {
            // For each diagonal element, the eigenvalues lie within radius equal to
            // the sum of absolute values of the adjacent off-diagonal elements
            let x = off_diagonal[i].abs() + off_diagonal[i + 1].abs();
            (diagonal[i] - x, diagonal[i] + x)
        })
        .reduce(
            // Initial value is the bounds for the last element
            || {
                (
                    diagonal[n - 1] - off_diagonal[n - 1].abs(),
                    diagonal[n - 1] + off_diagonal[n - 1].abs(),
                )
            },
            // Combine bounds by taking the minimum lower bound and maximum upper bound
            |(l1, u1), (l2, u2)| (min(l1, l2), max(u1, u2)),
        );

    (lower_bound, upper_bound)
}

// Computes the kth eigenvalue of a tridiagonal matrix using bisection
// k is 0-indexed, so k=0 means the smallest eigenvalue
#[time_graph::instrument]
fn kth_eigenvalue<F: Float + Send + Sync>(diagonal: &[F], off_diagonal: &[F], k: usize) -> F {
    let n = diagonal.len();
    let epsilon = F::from(f64::EPSILON).unwrap();
    let two = F::one() + F::one();

    // Get initial bounds for the eigenvalues
    let (mut xlower, mut xupper) = gershgorin_bounds(diagonal, off_diagonal);

    // Set initial tolerance based on the scale of the bounds
    let mut tolerance = two * epsilon * (xupper.abs() + xlower.abs());

    // Bisection method: repeatedly narrow the interval until desired precision is reached
    while (xupper - xlower).abs() > tolerance {
        let xmid = (xupper + xlower) / two;

        // Count how many eigenvalues are less than xmid
        let nb_eig_lt_xmid = nb_eigenvalues_lt_x(diagonal, off_diagonal, xmid);

        // Adjust bounds based on the count
        // We're looking for the kth eigenvalue from largest to smallest (n-k)
        if nb_eig_lt_xmid >= n - k {
            xupper = xmid;
        } else {
            xlower = xmid;
        }

        // Update tolerance based on current bounds
        tolerance = epsilon * (xupper.abs() + xlower.abs());
    }

    // Return the midpoint of the final interval
    (xlower + xupper) / two
}

// Computes all eigenvalues of a tridiagonal matrix in parallel
#[time_graph::instrument]
fn eigenvalues<F: Float + Send + Sync>(diagonal: &[F], off_diagonal: &[F]) -> Vec<F> {
    let n = diagonal.len();

    // Compute each eigenvalue in parallel
    let eigenvalues: Vec<F> = (0..n)
        .into_par_iter()
        .map(|k| kth_eigenvalue(diagonal, off_diagonal, k))
        .collect();

    eigenvalues
}

// Computes the zeros (roots) of the nth Laguerre polynomial
// Uses the fact that these zeros are eigenvalues of a specific tridiagonal matrix
#[time_graph::instrument]
pub fn laguerre_polynomial_zeros(n: usize) -> Vec<f64> {
    // Define the Jacobi matrix (tridiagonal symmetric matrix)
    // We first define the off-diagonal elements (these are √i in standard form)
    let off_diagonal: Vec<f64> = (0..n).into_par_iter().map(|i| i as f64).collect();

    // Then define the diagonal elements (these are 2i+1 in standard form)
    let diagonal: Vec<f64> = (0..n)
        .into_par_iter()
        .map(|i| {
            let e = 2 * i + 1;
            e as f64
        })
        .collect();

    // Compute the eigenvalues, which are the zeros of the Laguerre polynomial
    let zeros = eigenvalues(diagonal.as_slice(), off_diagonal.as_slice());
    return zeros;
}

// Test function to compute zeros of a large Laguerre polynomial
#[time_graph::instrument]
fn zeros() {
    let n = 1_000; // Degree of the polynomial
    laguerre_polynomial_zeros(n);
}

fn main() {
    // Enable performance profiling
    time_graph::enable_data_collection(true);

    // Compute the zeros
    zeros();

    // Get and print the performance profiling results
    let graph = time_graph::get_full_graph();
    // The following output formats are available but commented out:
    println!("{}", graph.as_dot()); // DOT format for visualization
                                    // println!("{}", graph.as_json());     // JSON format
                                    // println!("{}", graph.as_table());    // Full table
                                    // println!("{}", graph.as_short_table()); // Condensed table
}
