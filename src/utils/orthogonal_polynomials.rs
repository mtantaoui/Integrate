extern crate test;

use std::{fmt::Debug, ops::AddAssign};

use num::Float;

pub trait OrthogonalPolynomial<F: Float + Debug + AddAssign> {
    fn new(degree: usize) -> Self;

    fn eval(&self, x: F) -> F;

    fn eval_derivative(&self, x: F) -> F;

    fn zeros(&self) -> Vec<F>;
}
