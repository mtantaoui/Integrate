use std::fmt::Debug;
use std::{marker::PhantomData, ops::AddAssign};

use num::{Float, One, Zero};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::utils::matrix::TridiagonalSymmetricFloatMatrix;
use crate::utils::orthogonal_polynomials::OrthogonalPolynomial;

#[derive(Clone, Debug)]
pub struct Hermite<F: Float> {
    degree: usize,
    _x: PhantomData<F>,
}

impl<F: Float + Sync + Send + AddAssign + Debug> OrthogonalPolynomial<F> for Hermite<F> {
    fn new(degree: usize) -> Self {
        Hermite {
            degree,
            _x: PhantomData,
        }
    }

    fn eval(&self, x: F) -> F {
        let two = F::one() + F::one();
        let one = F::one();

        if self.degree.is_zero() {
            return F::one();
        }

        if self.degree.is_one() {
            return two * x;
        }

        let mut h_k_1 = F::one(); // H_{k-1}
        let mut h_k = two * x; // H_k

        let mut h = F::nan();

        for k in 2..=self.degree {
            let k = F::from(k).unwrap();

            h = two * x * h_k - two * (k - one) * h_k_1; // L_{k+1}

            h_k_1 = h_k;
            h_k = h;
        }

        h
    }

    fn zeros(&self) -> Vec<F> {
        let one = F::one();
        let two = F::one() + F::one();

        let diagonal = vec![F::zero(); self.degree];

        // let mut offdiagonal = vec![one / two; self.degree - 1];

        let mut offdiagonal: Vec<F> = (1..self.degree)
            .into_par_iter()
            .map(|i| {
                let i = F::from(i).unwrap();
                ((i + F::one()) / two).sqrt()
            })
            .collect();

        offdiagonal.insert(0, F::zero());
        let matrix = TridiagonalSymmetricFloatMatrix::new(diagonal, offdiagonal);

        matrix.eigenvalues()
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::FRAC_1_SQRT_2;

    use crate::{
        gauss_quadrature::hermite::Hermite, utils::orthogonal_polynomials::OrthogonalPolynomial,
    };

    const N_VALUES: &[usize; 18] = &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 5, 5, 5, 5, 5];
    const X_VALUES: &[f64; 18] = &[
        5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 0.0, 0.5, 1.0, 3.0, 10.0,
    ];
    const H_N_X: &[f64; 18] = &[
        1.0,
        10.0,
        98.0,
        940.0,
        8812.0,
        80600.0,
        717880.0,
        6211600.0,
        52065680.0,
        421271200.0,
        3275529760.0,
        24329873600.0,
        171237081280.0,
        0.0,
        41.0,
        -8.0,
        3816.0,
        3041200.0,
    ];

    const H1_ZEROS: [f64; 1] = [0.000000];
    // const H2_ZEROS: [f64; 2] = [-0.707107, 0.707107];
    const H2_ZEROS: [f64; 2] = [-FRAC_1_SQRT_2, FRAC_1_SQRT_2];
    const H3_ZEROS: [f64; 3] = [-1.224745, -0.000000, 1.224745];
    const H4_ZEROS: [f64; 4] = [-1.650680, -0.524648, 0.524648, 1.650680];
    const H5_ZEROS: [f64; 5] = [-2.020183, -0.958572, 0.000000, 0.958572, 2.020183];

    #[test]
    fn test_eval_laguerre() {
        for ((n, x), h) in N_VALUES.iter().zip(X_VALUES).zip(H_N_X) {
            let h_n = Hermite::new(*n);

            let h_n_x = h_n.eval(*x);

            assert_eq!(*h, h_n_x);
        }
    }

    #[test]
    fn test_hermite_zeros() {
        let h5: Hermite<f64> = Hermite::new(5);
        // let h5_zeros = h5.zeros();
        println!("{:?}", h5.zeros())
        // assert_eq!(h1_zeros, H1_ZEROS);
    }
}
