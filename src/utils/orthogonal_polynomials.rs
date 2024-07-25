extern crate test;

use core::f64;
use std::ops::AddAssign;

use num::{zero, Float, One, ToPrimitive, Zero};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::gauss_quadrature::bessel::bessel_j0_zeros;

pub trait OrthogonalPolynomial<F: Float> {
    fn new(degree: usize) -> Self;
    fn eval(&self, x: F) -> F;
    fn eval_derivative(&self, x: F) -> F;
    fn approximate_zero(&self, m: usize) -> f64;
    fn approximate_zeros(&self) -> Vec<F>;
    fn zeros(&self) -> Vec<F>;
    fn newton_raphson(&self, a: F, tolerance: F) -> F;
}

#[derive(Clone, Debug)]
pub struct LaguerrePolynomial<F: Float> {
    degree: usize,
    x: F,
}

impl<F: Float + Sync + Send + AddAssign> OrthogonalPolynomial<F> for LaguerrePolynomial<F> {
    fn new(degree: usize) -> LaguerrePolynomial<F> {
        LaguerrePolynomial {
            degree,
            x: F::nan(),
        }
    }

    fn eval(&self, x: F) -> F {
        if self.degree.is_zero() {
            return F::one();
        }

        if self.degree.is_one() {
            return F::one() - x;
        }

        let mut l_k_1 = F::one(); // L_{k-1}
        let mut l_k = F::one() - x; // L_k

        let mut l = F::nan();

        for k in 2..=self.degree {
            let a = F::from(2 * (k - 1) + 1).unwrap();
            let b = F::from(k - 1).unwrap();
            let c = F::from(k).unwrap();

            l = ((a - x) * l_k - b * l_k_1) / c; // L_{k+1}

            l_k_1 = l_k;
            l_k = l;
        }

        l
    }

    fn eval_derivative(&self, x: F) -> F {
        if x.is_zero() || self.degree.is_zero() {
            return zero();
        }

        let l_n_x = self.eval(x); // L_n(x)
        let l_n_1_x = LaguerrePolynomial::new(self.degree - 1).eval(x); // L_{n-1}(x)

        let n = F::from(self.degree).unwrap(); // converting n to Float to compute derivative

        (n * l_n_x - n * l_n_1_x) / x
    }

    fn approximate_zero(&self, m: usize) -> f64 {
        let n_f = self.degree.to_f64().unwrap();

        let j_0_m = bessel_j0_zeros(m);
        let k_n: f64 = n_f + 0.5;

        let term1 = j_0_m.powi(2) / (4.0 * k_n);
        let term2 = 1.0 + (j_0_m.powi(2) - 2.0) / (48.0 * k_n.powi(2));

        term1 * term2
    }

    fn approximate_zeros(&self) -> Vec<F> {
        (1..=self.degree)
            .into_par_iter()
            .map(|m| {
                let a = self.approximate_zero(m);
                F::from(a).unwrap()
            })
            .collect()
    }

    fn zeros(&self) -> Vec<F> {
        let approx = self.approximate_zeros();

        let tolerance = F::from(10e-3).unwrap();

        approx
            .into_par_iter()
            .map(|a| self.newton_raphson(a, tolerance))
            .collect()
    }

    fn newton_raphson(&self, mut a: F, tolerance: F) -> F {
        let mut delta;
        let mut num_iterations: usize = 0;

        loop {
            let fa = self.eval(a);

            if fa.is_zero() {
                return a;
            }

            let dfa = self.eval_derivative(a);

            // If an attempt to divide by zero, TODO: return a specific error
            if dfa.is_zero() {
                panic!()
            }

            delta = -fa / dfa;

            // println!("{}", delta.to_f64().unwrap());

            // If the location of the root relative to a is less than //
            // the tolerance, return the estimate of the root.        //

            if delta.abs() < tolerance || num_iterations > 100 {
                break;
            }

            a += delta;
            num_iterations += 1;
        }

        println!("{}", num_iterations);
        a + delta
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use test::Bencher;

    #[bench]
    fn bench_laguerre_zeros(bencher: &mut Bencher) {
        let lag: LaguerrePolynomial<f64> = LaguerrePolynomial::new(1_000);
        bencher.iter(|| lag.zeros())
    }
}
