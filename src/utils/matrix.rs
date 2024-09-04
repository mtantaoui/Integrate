extern crate test;

use num::Float;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

pub struct TridiagonalSymmetricFloatMatrix<F: Float> {
    diagonal: Vec<F>,
    offdiagonal: Vec<F>,
}

impl<F: Float + Send + Sync> TridiagonalSymmetricFloatMatrix<F> {
    pub fn new(diagonal: Vec<F>, offdiagonal: Vec<F>) -> TridiagonalSymmetricFloatMatrix<F> {
        TridiagonalSymmetricFloatMatrix {
            diagonal,
            offdiagonal,
        }
    }

    pub fn eigenvalues(&self) -> Vec<F> {
        let n = self.diagonal.len();
        let eigenvalues: Vec<F> = (0..n)
            .into_par_iter()
            .map(|k| self.kth_eigenvalue(k))
            .collect();
        eigenvalues
    }

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

            if tolerance.is_zero() {
                tolerance = epsilon
            }
        }

        (xlower + xupper) / two
    }

    fn gershgorin_bounds(&self) -> (F, F) {
        let n = self.diagonal.len();

        let (lower_bound, upper_bound) = (0..n - 1)
            .into_par_iter()
            .map(|i| {
                let x = self.offdiagonal[i].abs() + self.offdiagonal[i + 1].abs();
                (self.diagonal[i] - x, self.diagonal[i] + x)
            })
            .reduce(
                || {
                    (
                        self.diagonal[n - 1] - self.offdiagonal[n - 1].abs(),
                        self.diagonal[n - 1] + self.offdiagonal[n - 1].abs(),
                    )
                },
                |(l1, u1), (l2, u2)| (min(l1, l2), max(u1, u2)),
            );

        (lower_bound, upper_bound)
    }

    fn nb_eigenvalues_lt_x(&self, x: F) -> usize {
        let mut q = F::one();
        let epsilon = F::from(f64::EPSILON).unwrap();
        let mut k: usize = 0;
        let n = self.diagonal.len();

        for i in 0..n {
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

fn max<F: Float>(a: F, b: F) -> F {
    let two: F = F::one() + F::one();
    ((a + b) + (a - b).abs()) / two
}

fn min<F: Float>(a: F, b: F) -> F {
    let two: F = F::one() + F::one();
    ((a + b) - (a - b).abs()) / two
}

#[cfg(test)]
mod tests {
    use super::*;
    use test::Bencher;

    #[test]
    fn test_tdsf_matrix() {
        let n: usize = 1_000;
        let diagonal: Vec<f64> = (1..=n).map(|e| e.pow(2) as f64).collect();
        let offdiagonal: Vec<f64> = (0..n).map(|e| e.pow(4) as f64).collect();
        let matrix = TridiagonalSymmetricFloatMatrix::new(diagonal, offdiagonal);

        matrix.eigenvalues();
    }

    #[bench]
    fn bench_eigenvalues(bencher: &mut Bencher) {
        let n: usize = 1_000;
        let diagonal: Vec<f64> = (1..=n).map(|e| e.pow(2) as f64).collect();
        let offdiagonal: Vec<f64> = (0..n).map(|e| e.pow(4) as f64).collect();
        let matrix = TridiagonalSymmetricFloatMatrix::new(diagonal, offdiagonal);

        bencher.iter(|| {
            matrix.eigenvalues();
        })
    }
}
