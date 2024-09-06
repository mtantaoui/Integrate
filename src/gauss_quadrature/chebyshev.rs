use std::iter::Sum;
use std::{f64::consts::PI, marker::PhantomData};

use std::fmt::Debug;
use std::ops::AddAssign;

use num::{one, Float, Zero};
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};

use crate::utils::orthogonal_polynomials::OrthogonalPolynomial;

use super::utils::check_gauss_rule_args;

#[derive(Clone, Debug)]
pub struct ChebyshevFirstKind<F: Float> {
    degree: usize,
    _x: PhantomData<F>,
}

#[derive(Clone, Debug)]
pub struct ChebyshevSecondKind<F: Float> {
    degree: usize,
    _x: PhantomData<F>,
}

fn roots_first_kind_chebyshev<F: Float + Debug + Sync + Send + AddAssign>(
    n: usize,
) -> (Vec<F>, Vec<F>) {
    let t_n: ChebyshevFirstKind<F> = ChebyshevFirstKind::new(n);
    let zeros = t_n.zeros();

    let pi = F::from(PI).unwrap();
    let n = F::from(n).unwrap();

    let weights = vec![pi / n; t_n.degree];

    (zeros, weights)
}

fn roots_second_kind_chebyshev<F: Float + Debug + Sync + Send + AddAssign>(
    n: usize,
) -> (Vec<F>, Vec<F>) {
    let u_n: ChebyshevSecondKind<F> = ChebyshevSecondKind::new(n);
    let zeros = u_n.zeros();

    let n = F::from(u_n.degree).unwrap();
    let pi = F::from(PI).unwrap();

    let weights: Vec<F> = (1..=u_n.degree)
        .into_par_iter()
        .map(|i| {
            let i = F::from(i).unwrap();

            let numer = i * pi;
            let denom = n + one();

            let angle = numer / denom;

            let term1 = angle.sin().powi(2);
            let term2 = pi / (n + one());

            term1 * term2
        })
        .collect();

    (zeros, weights)
}

pub fn gauss_first_kind_chebyshev_rule<F: Float + Debug + Sync + Send + AddAssign + Sum>(
    f: fn(F) -> F,
    n: usize,
) -> F {
    check_gauss_rule_args(n);

    let (zeros, weights) = roots_first_kind_chebyshev::<F>(n);

    weights
        .into_par_iter()
        .zip(zeros)
        .map(|(w, x)| w * f(x))
        .sum()
}

pub fn gauss_second_kind_chebyshev_rule<F: Float + Debug + Sync + Send + AddAssign + Sum>(
    f: fn(F) -> F,
    n: usize,
) -> F {
    check_gauss_rule_args(n);

    let (zeros, weights) = roots_second_kind_chebyshev::<F>(n);

    weights
        .into_par_iter()
        .zip(zeros)
        .map(|(w, x)| w * f(x))
        .sum()
}

impl<F: Float + Debug + AddAssign + Send + Sync> OrthogonalPolynomial<F> for ChebyshevFirstKind<F> {
    fn new(degree: usize) -> Self {
        ChebyshevFirstKind {
            degree,
            _x: PhantomData,
        }
    }

    fn eval(&self, x: F) -> F {
        let theta = x.acos();
        let n = F::from(self.degree).unwrap();

        (n * theta).cos()
    }

    fn zeros(&self) -> Vec<F> {
        if self.degree.is_zero() {
            return vec![];
        }

        let n = F::from(self.degree).unwrap();
        let pi = F::from(PI).unwrap();
        let two = F::one() + F::one();

        let zeros: Vec<F> = (1..=self.degree)
            .into_par_iter()
            .map(|i| {
                let i = F::from(i).unwrap();

                let numer = (two * i - one()) * pi;
                let denom = two * n;

                let angle = numer / denom;

                angle.cos()
            })
            .collect();

        zeros
    }
}

impl<F: Float + Debug + AddAssign + Send + Sync> OrthogonalPolynomial<F>
    for ChebyshevSecondKind<F>
{
    fn new(degree: usize) -> Self {
        ChebyshevSecondKind {
            degree,
            _x: PhantomData,
        }
    }

    fn eval(&self, x: F) -> F {
        let theta = x.acos();
        let n = F::from(self.degree).unwrap();

        let numer = ((n + one()) * theta).sin();
        let denom = theta.sin();

        numer / denom
    }

    fn zeros(&self) -> Vec<F> {
        if self.degree.is_zero() {
            return vec![];
        }

        let n = F::from(self.degree).unwrap();
        let pi = F::from(PI).unwrap();

        let zeros: Vec<F> = (1..=self.degree)
            .into_par_iter()
            .map(|i| {
                let i = F::from(i).unwrap();

                let numer = i * pi;
                let denom = n + one();

                let angle = numer / denom;

                angle.cos()
            })
            .collect();

        zeros
    }
}
