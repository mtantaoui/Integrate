use std::{f64::consts::PI, marker::PhantomData};

use std::fmt::Debug;
use std::ops::AddAssign;

use num::{one, Float, Zero};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::utils::orthogonal_polynomials::OrthogonalPolynomial;

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
