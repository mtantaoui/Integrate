use std::{
    f64::consts::PI,
    fmt::Debug,
    marker::PhantomData,
    ops::{AddAssign, Mul},
};

use num::{bigint::ToBigInt, BigRational, BigUint, Float, One, Zero};
use rayon::iter::{
    IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator, ParallelExtend,
    ParallelIterator,
};

use super::{
    matrix::TridiagonalSymmetricFloatMatrix, orthogonal_polynomials::OrthogonalPolynomial,
};

#[derive(Clone, Debug)]
struct Hermite<F: Float> {
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
        if self.degree.is_zero() {
            return vec![];
        }

        let two = F::one() + F::one();

        // define the Jacobi matrix (tridiagonal symmetric matrix)
        let diagonal = vec![F::zero(); self.degree];

        let mut offdiagonal = vec![F::zero()];
        offdiagonal.par_extend((0..self.degree - 1).into_par_iter().map(|i| {
            let i = F::from(i).unwrap();
            ((i + F::one()) / two).sqrt()
        }));

        let matrix = TridiagonalSymmetricFloatMatrix::new(diagonal, offdiagonal);

        matrix.eigenvalues()
    }
}

// weights formula : https://wikimedia.org/api/rest_v1/media/math/render/svg/2e6f152a1e9ecd4ab8ddf912aaa69bb8d0e66a3c
pub fn roots_hermite<F: Float + Debug + AddAssign + Sync + Send + ToBigInt>(
    n: usize,
) -> (Vec<F>, Vec<F>) {
    let h_n: Hermite<F> = Hermite::new(n); // H_n
    let zeros = h_n.zeros();

    let h: Hermite<F> = Hermite::new(n - 1); // H_{n-1}

    // params used in weights formula
    let sqrt_pi = F::from(PI).unwrap().sqrt();

    let n_fact = F::from(factorial(n)).unwrap();

    let two = F::one() + F::one();
    let n_squared = F::from(n).unwrap().powf(two);
    let n = F::from(n).unwrap();

    let two_pow = two.powf(n - F::one());

    let weights: Vec<F> = zeros
        .par_iter()
        .map(|x_i| {
            let h_x = h.eval(*x_i); // H_{n-1}(x_i)

            let numerator = two_pow * n_fact * sqrt_pi;

            let denominator = n_squared * h_x * h_x;

            if denominator.is_infinite() || numerator.is_infinite() {
                // switching everything to BigInt
                let numer = two_pow.to_bigint().unwrap() * n_fact.to_bigint().unwrap();
                let denom = h_x.abs().to_bigint().unwrap().pow(2) * n_squared.to_bigint().unwrap();
                let ratio = BigRational::new(numer, denom);

                F::from(ratio).unwrap() * sqrt_pi
            } else {
                numerator / denominator
            }
        })
        .collect();

    let warn = zeros
        .as_slice()
        .into_par_iter()
        .zip(weights.as_slice())
        .any(|(zero, weight)| (*zero).is_nan() || (*weight).is_nan());

    if warn {
        eprintln!(
            "Warning: `n` chosen is too big, some values of Hermite Polynomials weights or zeros are too small and may underflow!"
        )
    }

    (zeros, weights)
}

fn factorial(n: usize) -> BigUint {
    (1..n + 1)
        .into_par_iter()
        // .with_min_len(64)
        .fold_with(BigUint::from(1_usize), |acc, x| acc.mul(x))
        .reduce_with(Mul::mul)
        .unwrap()
}

#[cfg(test)]
mod tests {
    use std::f64::consts::FRAC_1_SQRT_2;

    use crate::utils::{hermite::Hermite, orthogonal_polynomials::OrthogonalPolynomial};

    const EPSILON: f64 = 10e-7;

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
        let h1: Hermite<f64> = Hermite::new(1);
        let h2: Hermite<f64> = Hermite::new(2);
        let h3: Hermite<f64> = Hermite::new(3);
        let h4: Hermite<f64> = Hermite::new(4);
        let h5: Hermite<f64> = Hermite::new(5);

        let mut h1_zeros = h1.zeros();
        h1_zeros.reverse();

        let mut h2_zeros = h2.zeros();
        h2_zeros.reverse();

        let mut h3_zeros = h3.zeros();
        h3_zeros.reverse();

        let mut h4_zeros = h4.zeros();
        h4_zeros.reverse();

        let mut h5_zeros = h5.zeros();
        h5_zeros.reverse();

        let h1_test = h1_zeros
            .iter()
            .zip(H1_ZEROS)
            .all(|(value, test_value)| (value - test_value).abs() <= EPSILON);

        assert!(h1_test);

        let h2_test = h2_zeros
            .iter()
            .zip(H2_ZEROS)
            .all(|(value, test_value)| (value - test_value).abs() <= EPSILON);

        assert!(h2_test);

        let h3_test = h3_zeros
            .iter()
            .zip(H3_ZEROS)
            .all(|(value, test_value)| (value - test_value).abs() <= EPSILON);

        assert!(h3_test);

        let h4_test = h4_zeros
            .iter()
            .zip(H4_ZEROS)
            .all(|(value, test_value)| (value - test_value).abs() <= EPSILON);

        assert!(h4_test);

        let h5_test = h5_zeros
            .iter()
            .zip(H5_ZEROS)
            .all(|(value, test_value)| (value - test_value).abs() <= EPSILON);

        assert!(h5_test)
    }
}
