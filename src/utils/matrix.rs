// extern crate test;

use num::Float;

/// A symmetric tridiagonal matrix of floating-point values.
///
/// Used by the Golub-Welsch algorithm to compute the zeros of orthogonal
/// polynomials: the zeros coincide with the eigenvalues of the Jacobi
/// tridiagonal matrix built from the polynomial's recurrence coefficients.
pub struct TridiagonalSymmetricFloatMatrix<F: Float> {
    diagonal: Vec<F>,
    offdiagonal: Vec<F>,
}

impl<F: Float + Send + Sync> TridiagonalSymmetricFloatMatrix<F> {
    /// Creates a new symmetric tridiagonal matrix.
    ///
    /// * `diagonal`    - the $n$ main-diagonal entries $\alpha_1, \ldots, \alpha_n$.
    /// * `offdiagonal` - the $n$ sub/super-diagonal entries; `offdiagonal[0]` is
    ///   unused (treated as zero) and `offdiagonal[i]` holds $\beta_i$ for $i \ge 1$.
    pub fn new(diagonal: Vec<F>, offdiagonal: Vec<F>) -> TridiagonalSymmetricFloatMatrix<F> {
        TridiagonalSymmetricFloatMatrix {
            diagonal,
            offdiagonal,
        }
    }

    /// Returns all eigenvalues of the matrix.
    ///
    /// Each eigenvalue is found independently via [`Self::kth_eigenvalue`].
    pub fn eigenvalues(&self) -> Vec<F> {
        let n = self.diagonal.len();
        let eigenvalues: Vec<F> = (0..n).into_iter().map(|k| self.kth_eigenvalue(k)).collect();
        eigenvalues
    }

    /// Returns the $k$-th smallest eigenvalue using bisection on Sturm sequences.
    ///
    /// The algorithm maintains a bracket `[xlower, xupper]` initialized from
    /// [`Self::gershgorin_bounds`] and halves it each iteration: after
    /// computing the midpoint, [`Self::nb_eigenvalues_lt_x`] counts how many
    /// eigenvalues lie below the midpoint. If that count is $\ge n - k$ the
    /// upper bound is moved down; otherwise the lower bound is moved up.
    /// The loop terminates when the bracket width falls within floating-point
    /// tolerance relative to the current bounds.
    fn kth_eigenvalue(&self, k: usize) -> F {
        let n = self.diagonal.len();
        let epsilon = F::from(f32::EPSILON).unwrap();
        let two = F::one() + F::one();

        let (mut xlower, mut xupper) = self.gershgorin_bounds();

        let mut tolerance = two * epsilon * (xupper.abs() + xlower.abs());

        while (xupper - xlower).abs() > tolerance {
            let xmid = (xupper + xlower) / two;

            let nb_eig_lt_xmid = self.nb_eigenvalues_lt_x(xmid);

            if nb_eig_lt_xmid >= n - k {
                xupper = xmid;
            } else {
                xlower = xmid;
            }

            tolerance = epsilon * (xupper.abs() + xlower.abs());

            // Avoid an infinite loop when the bracket collapses to zero.
            if tolerance.is_zero() {
                tolerance = epsilon
            }
        }

        (xlower + xupper) / two
    }

    /// Returns a `(lower, upper)` interval guaranteed to contain all eigenvalues,
    /// derived from the Gershgorin circle theorem.
    ///
    /// For each row $i$ the Gershgorin disk has center $\alpha_i$ and radius
    /// equal to the sum of the absolute off-diagonal entries in that row.
    /// The returned bounds are the tightest such interval over all rows.
    fn gershgorin_bounds(&self) -> (F, F) {
        let n = self.diagonal.len();

        // Start with the Gershgorin disk for the last row (only one off-diagonal entry).
        let init = (
            self.diagonal[n - 1] - self.offdiagonal[n - 1].abs(),
            self.diagonal[n - 1] + self.offdiagonal[n - 1].abs(),
        );

        let (lower_bound, upper_bound) = (0..n - 1)
            .map(|i| {
                // Radius of the Gershgorin disk for row i.
                let x = self.offdiagonal[i].abs() + self.offdiagonal[i + 1].abs();
                (self.diagonal[i] - x, self.diagonal[i] + x)
            })
            .fold(init, |(l1, u1), (l2, u2)| (min(l1, l2), max(u1, u2)));

        (lower_bound, upper_bound)
    }

    /// Returns the number of eigenvalues strictly less than `x`.
    ///
    /// Uses the Sturm sequence recurrence: starting from $q_0 = 1$, each term
    /// is updated as
    ///
    /// $$q_i = (\alpha_i - x) - \frac{\beta_i^2}{q_{i-1}}$$
    ///
    /// where a zero denominator is replaced by a small $\varepsilon$ to avoid
    /// division by zero. The count of sign changes in the sequence equals the
    /// number of eigenvalues less than `x`.
    fn nb_eigenvalues_lt_x(&self, x: F) -> usize {
        let mut q = F::one();
        let epsilon = F::from(f64::EPSILON).unwrap();
        let mut k: usize = 0;
        let n = self.diagonal.len();

        for i in 0..n {
            // Sturm sequence recurrence; guard against q == 0.
            q = if q.is_zero() {
                self.diagonal[i] - x - self.offdiagonal[i].abs() / epsilon
            } else {
                self.diagonal[i] - x - self.offdiagonal[i].powi(2) / q
            };

            if q.is_sign_negative() {
                k += 1;
            }
        }

        k
    }
}

/// Returns the larger of `a` and `b` without branching.
fn max<F: Float>(a: F, b: F) -> F {
    let two: F = F::one() + F::one();
    ((a + b) + (a - b).abs()) / two
}

/// Returns the smaller of `a` and `b` without branching.
fn min<F: Float>(a: F, b: F) -> F {
    let two: F = F::one() + F::one();
    ((a + b) - (a - b).abs()) / two
}

#[cfg(test)]
mod tests {
    use super::*;
    // use test::Bencher;

    #[test]
    fn test_tdsf_matrix() {
        let n: usize = 1_000;
        let diagonal: Vec<f64> = (1..=n).map(|e| e.pow(2) as f64).collect();
        let offdiagonal: Vec<f64> = (0..n).map(|e| e.pow(4) as f64).collect();
        let matrix = TridiagonalSymmetricFloatMatrix::new(diagonal, offdiagonal);

        matrix.eigenvalues();
    }

    // #[bench]
    // fn bench_eigenvalues(bencher: &mut Bencher) {
    //     let n: usize = 1_000;
    //     let diagonal: Vec<f64> = (1..=n).map(|e| e.pow(2) as f64).collect();
    //     let offdiagonal: Vec<f64> = (0..n).map(|e| e.pow(4) as f64).collect();
    //     let matrix = TridiagonalSymmetricFloatMatrix::new(diagonal, offdiagonal);

    //     bencher.iter(|| {
    //         matrix.eigenvalues();
    //     })
    // }
}
