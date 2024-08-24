extern crate test;

use std::{fmt::Debug, marker::PhantomData, ops::AddAssign};

use num::{one, zero, Float, One, ToPrimitive, Zero};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::gauss_quadrature::bessel::bessel_j0_zeros;

const MAX_ITERATIONS: usize = 1_000;

pub trait OrthogonalPolynomial<F: Float + Debug + AddAssign> {
    fn new(degree: usize) -> Self;

    fn eval(&self, x: F) -> F;

    fn eval_derivative(&self, x: F) -> F;

    fn approximate_zero(&self, m: usize) -> f64;

    fn approximate_zeros(&self) -> Vec<F>;

    fn zero_upper_bound(&self, m: usize) -> F;

    fn zero_lower_bound(&self, m: usize) -> F;

    fn zeros_upper_bounds(&self) -> Vec<F>;

    fn zeros_lower_bounds(&self) -> Vec<F>;

    fn zeros(&self) -> Vec<F>;

    fn zeros_amsterdam(&self) -> Vec<F>;

    fn newton_raphson(&self, a: F, tolerance: F) -> F {
        let mut a = a;
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

            if delta.abs() < tolerance || num_iterations > MAX_ITERATIONS {
                break;
            }

            a += delta;
            num_iterations += 1;
        }

        a
    }

    fn amsterdam(&self, a: F, c: F, tolerance: F) -> F {
        let two = F::one() + one();

        let mut a = a;
        let mut c = c;

        let mut fa = self.eval(a);
        let mut fc = self.eval(c);
        let mut fb = fa * fc;

        let mut b = (a + c) / two;

        let mut delta: F;
        let mut dab: F;
        let mut dcb: F;

        // If the initial estimates do not bracket a root, set the err flag and
        // return.  If an initial estimate is a root, then return the root.

        if fb >= zero() {
            if fb > zero() {
                let err_msg = format!(
                    "root outside the passed interval [a, c] = [{:?}, {:?}]\nf(a)={:?}\nf(c)={:?}",
                    a, c, fa, fc
                );
                panic!("{}", err_msg);
            } else {
                let root = if fa.is_zero() { a } else { c };
                return root;
            }
        }

        // Insure that the initial estimate a < c.

        if a > c {
            delta = a;
            a = c;
            c = delta;
            delta = fa;
            fa = fc;
            fc = delta;
        }

        // If the function at the left endpoint is positive, and the function
        // at the right endpoint is negative.  Iterate reducing the length
        // of the interval by either bisection or quadratic inverse
        // interpolation.  Note that the function at the left endpoint
        // remains nonnegative and the function at the right endpoint remains
        // nonpositive.

        if fa > zero() {
            for _ in 0..MAX_ITERATIONS {
                // Are the two endpoints within the user specified tolerance ?

                if c - a < tolerance {
                    return (a + c) / two;
                };

                // No, Continue iteration.

                fb = self.eval(b);

                // Check that we are converging or that we have converged near
                // the left endpoint of the inverval.  If it appears that the
                // interval is not decreasing fast enough, use bisection.

                if b - a < tolerance {
                    if fb > zero() {
                        a = b;
                        fa = fb;
                        b = (a + c) / two;
                        continue;
                    } else {
                        return b;
                    }
                }

                // Check that we are converging or that we have converged near
                // the right endpoint of the inverval.  If it appears that the
                // interval is not decreasing fast enough, use bisection.

                if c - b < tolerance {
                    if fb < zero() {
                        c = b;
                        fc = fb;
                        b = (a + c) / two;
                        continue;
                    } else {
                        return b;
                    }
                }

                // If quadratic inverse interpolation is feasible, try it.
                if fa > fb && fb > fc {
                    delta = denominator(fa, fb, fc);

                    if delta != zero() {
                        dab = a - b;
                        dcb = c - b;
                        delta = numerator(dab, dcb, fa, fb, fc) / delta;

                        // Will the new estimate of the root be within the
                        // interval?  If yes, use it and update interval.
                        // If no, use the bisection method.

                        if delta > dab && delta < dcb {
                            if fb > zero() {
                                a = b;
                                fa = fb;
                            } else if fb < zero() {
                                c = b;
                                fc = fb;
                            } else {
                                return b;
                            }
                            b += delta;
                            continue;
                        }
                    }
                }

                if fb > zero() {
                    a = b;
                    fa = fb;
                } else {
                    c = b;
                    fc = fb;
                }

                b = (a + c) / two;
            }
        } else {
            // If the function at the left endpoint is negative, and the function
            // at the right endpoint is positive.  Iterate reducing the length
            // of the interval by either bisection or quadratic inverse
            // interpolation.  Note that the function at the left endpoint
            // remains nonpositive and the function at the right endpoint remains
            // nonnegative.

            for _ in 0..MAX_ITERATIONS {
                if c - a < tolerance {
                    return (a + c) / two;
                };

                fb = self.eval(b);

                if b - a < tolerance {
                    if fb < zero() {
                        a = b;
                        fa = fb;
                        b = (a + c) / two;
                        continue;
                    } else {
                        return b;
                    }
                }

                if c - b < tolerance {
                    if fb > zero() {
                        c = b;
                        fc = fb;
                        b = (a + c) / two;
                        continue;
                    } else {
                        return b;
                    }
                }

                if fa < fb && fb < fc {
                    delta = denominator(fa, fb, fc);

                    if delta != zero() {
                        dab = a - b;
                        dcb = c - b;
                        delta = numerator(dab, dcb, fa, fb, fc) / delta;

                        if delta > dab && delta < dcb {
                            if fb < zero() {
                                a = b;
                                fa = fb;
                            } else if fb > zero() {
                                c = b;
                                fc = fb;
                            } else {
                                return b;
                            }
                            b += delta;
                            continue;
                        }
                    }
                }

                if fb < zero() {
                    a = b;
                    fa = fb;
                } else {
                    c = b;
                    fc = fb;
                }

                b = (a + c) / two;
            }
        }
        b
    }
}

#[derive(Clone, Debug)]
pub struct LaguerrePolynomial<F: Float> {
    degree: usize,
    _x: PhantomData<F>,
}

impl<F: Float + Sync + Send + AddAssign + Debug> OrthogonalPolynomial<F> for LaguerrePolynomial<F> {
    fn new(degree: usize) -> LaguerrePolynomial<F> {
        LaguerrePolynomial {
            degree,
            _x: PhantomData,
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
        // let k_n: f64 = n_f + 0.5;

        // let term1 = j_0_m.powi(2) / (4.0 * k_n);
        // let term2 = 1.0 + (j_0_m.powi(2) - 2.0) / (48.0 * k_n.powi(2));

        // term1 * term2
        j_0_m.powi(2) / (4.0 * n_f + 2.0)
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

        let tolerance = F::from(10e-6).unwrap();

        approx
            .into_par_iter()
            .map(|a| self.newton_raphson(a, tolerance))
            .collect()
    }

    fn zero_upper_bound(&self, m: usize) -> F {
        let quarter = F::from(0.25).unwrap();

        fn term<F: Float>(m: usize) -> F {
            let t = 2 * m + 1;
            F::from(t).unwrap()
        }

        let numerator: F = term::<F>(m) * (term::<F>(m) + (term::<F>(m).powi(2) + quarter).sqrt());
        let denominator: F = term::<F>(self.degree);

        numerator / denominator
    }

    fn zero_lower_bound(&self, m: usize) -> F {
        let bessel_j0_squared = bessel_j0_zeros(m).powi(2);

        let numerator = F::from(bessel_j0_squared).unwrap();
        let denominator: F = mu::<F>(self.degree);

        numerator / denominator
    }

    fn zeros_upper_bounds(&self) -> Vec<F> {
        (1..=self.degree)
            .into_par_iter()
            .map(|m| self.zero_upper_bound(m))
            .collect()
    }

    fn zeros_lower_bounds(&self) -> Vec<F> {
        (1..=self.degree)
            .into_par_iter()
            .map(|m| self.zero_lower_bound(m))
            .collect()
    }

    fn zeros_amsterdam(&self) -> Vec<F> {
        todo!()

        // let tolerance = F::from(10e-6).unwrap();

        // let lower_bounds: Vec<F> = (1..=self.degree)
        //     .into_par_iter()
        //     .map(|m| self.zero_lower_bound(m))
        //     .collect();

        // let upper_bounds: Vec<F> = (1..=self.degree)
        //     .into_par_iter()
        //     .map(|m| self.zero_upper_bound(m))
        //     .collect();

        // lower_bounds
        //     .into_par_iter()
        //     .zip(upper_bounds)
        //     .map(|(lower_bound, upper_bound)| self.amsterdam(lower_bound, upper_bound, tolerance))
        //     .collect()
    }
}

fn mu<F: Float>(n: usize) -> F {
    let mu = 4 * n + 2;
    F::from(mu).unwrap()
}

fn numerator<F: Float>(dab: F, dcb: F, fa: F, fb: F, fc: F) -> F {
    fb * (dab * fc * (fc - fb) - fa * dcb * (fa - fb))
}

fn denominator<F: Float>(fa: F, fb: F, fc: F) -> F {
    (fc - fb) * (fa - fb) * (fa - fc)
}

#[cfg(test)]
mod tests {

    use super::*;
    use test::Bencher;

    const EPSILON: f64 = 10e-10;

    const L_N_X: &[f64; 17] = &[
        0.1000000000000000E+01,
        0.0000000000000000E+00,
        -0.5000000000000000E+00,
        -0.6666666666666667E+00,
        -0.6250000000000000E+00,
        -0.4666666666666667E+00,
        -0.2569444444444444E+00,
        -0.4047619047619048E-01,
        0.1539930555555556E+00,
        0.3097442680776014E+00,
        0.4189459325396825E+00,
        0.4801341790925124E+00,
        0.4962122235082305E+00,
        -0.4455729166666667E+00,
        0.8500000000000000E+00,
        -0.3166666666666667E+01,
        0.3433333333333333E+02,
    ];

    const L_N_X_DERIV: &[f64; 17] = &[
        0.0,
        -1.0,
        -1.0,
        -0.5,
        0.16666666666666785,
        0.7916666666666679,
        1.258333333333331,
        1.5152777777777793,
        1.5557539682539598,
        1.4017609126984052,
        1.092016644620804,
        0.673070712081147,
        0.19293653298857372,
        -1.1484374999999991,
        -0.8749999999999991,
        -1.8749999999999991,
        11.666666666666785,
    ];

    const N_VALUES: &[usize; 17] = &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 5, 5, 5, 5];

    const X_VALUES: &[f64; 17] = &[
        1.0E+00, 1.0E+00, 1.0E+00, 1.0E+00, 1.0E+00, 1.0E+00, 1.0E+00, 1.0E+00, 1.0E+00, 1.0E+00,
        1.0E+00, 1.0E+00, 1.0E+00, 0.5E+00, 3.0E+00, 5.0E+00, 1.0E+01,
    ];

    // computed using scipy roots_laguerre
    const FIRST_100_LAGUERRE_ROOTS: [f64; 100] = [
        1.43861470e-02,
        7.58036120e-02,
        1.86314102e-01,
        3.45969181e-01,
        5.54810938e-01,
        8.12891284e-01,
        1.12027384e+00,
        1.47703433e+00,
        1.88326083e+00,
        2.33905385e+00,
        2.84452654e+00,
        3.39980483e+00,
        4.00502758e+00,
        4.66034684e+00,
        5.36592799e+00,
        6.12195003e+00,
        6.92860583e+00,
        7.78610238e+00,
        8.69466111e+00,
        9.65451824e+00,
        1.06659251e+01,
        1.17291485e+01,
        1.28444712e+01,
        1.40121922e+01,
        1.52326276e+01,
        1.65061105e+01,
        1.78329919e+01,
        1.92136416e+01,
        2.06484480e+01,
        2.21378194e+01,
        2.36821848e+01,
        2.52819939e+01,
        2.69377187e+01,
        2.86498542e+01,
        3.04189188e+01,
        3.22454560e+01,
        3.41300351e+01,
        3.60732524e+01,
        3.80757324e+01,
        4.01381291e+01,
        4.22611275e+01,
        4.44454451e+01,
        4.66918335e+01,
        4.90010802e+01,
        5.13740103e+01,
        5.38114889e+01,
        5.63144230e+01,
        5.88837640e+01,
        6.15205103e+01,
        6.42257101e+01,
        6.70004645e+01,
        6.98459306e+01,
        7.27633254e+01,
        7.57539295e+01,
        7.88190913e+01,
        8.19602322e+01,
        8.51788512e+01,
        8.84765308e+01,
        9.18549433e+01,
        9.53158573e+01,
        9.88611460e+01,
        1.02492795e+02,
        1.06212911e+02,
        1.10023736e+02,
        1.13927651e+02,
        1.17927199e+02,
        1.22025092e+02,
        1.26224231e+02,
        1.30527723e+02,
        1.34938905e+02,
        1.39461365e+02,
        1.44098970e+02,
        1.48855901e+02,
        1.53736688e+02,
        1.58746249e+02,
        1.63889946e+02,
        1.69173640e+02,
        1.74603761e+02,
        1.80187391e+02,
        1.85932360e+02,
        1.91847369e+02,
        1.97942133e+02,
        2.04227560e+02,
        2.10715973e+02,
        2.17421393e+02,
        2.24359895e+02,
        2.31550068e+02,
        2.39013630e+02,
        2.46776241e+02,
        2.54868629e+02,
        2.63328168e+02,
        2.72201170e+02,
        2.81546328e+02,
        2.91440134e+02,
        3.01985855e+02,
        3.13329534e+02,
        3.25691263e+02,
        3.39435102e+02,
        3.55261312e+02,
        3.74984113e+02,
    ];

    #[test]
    fn test_eval_laguerre() {
        for ((&ln_test, &n), &x) in L_N_X.into_iter().zip(N_VALUES).zip(X_VALUES) {
            let lag: LaguerrePolynomial<f64> = LaguerrePolynomial::new(n);

            let ln = lag.eval(x);

            let is_close = (ln - ln_test).abs() < EPSILON;

            assert!(is_close)
        }
    }

    #[test]
    fn test_eval_laguerre_derivative() {
        for ((&dln_test, &n), &x) in L_N_X_DERIV.into_iter().zip(N_VALUES).zip(X_VALUES) {
            let lag: LaguerrePolynomial<f64> = LaguerrePolynomial::new(n);

            let dln = lag.eval_derivative(x);

            let is_close = (dln - dln_test).abs() < EPSILON;

            assert!(is_close)
        }
    }

    #[test]
    fn test_laguerre_polynomial_zeros() {
        todo!()
        // let lag: LaguerrePolynomial<f64> = LaguerrePolynomial::new(100);

        // let lag_zeros = lag.zeros();

        // for (test_zero, zero) in lag_zeros.into_iter().zip(FIRST_100_LAGUERRE_ROOTS) {
        //     // println!("test: {}\t zero:{}", lag.eval(test_zero), lag.eval(zero));
        //     println!("difference:{}", lag.eval(test_zero) < lag.eval(zero));
        // }

        // // assert!(is_close)
        // // assert!(true);
    }

    #[bench]
    fn bench_laguerre_zeros_newton_raphson(bencher: &mut Bencher) {
        let lag: LaguerrePolynomial<f64> = LaguerrePolynomial::new(1_000);
        bencher.iter(|| lag.zeros())
    }

    #[bench]
    fn bench_laguerre_zeros_amsterdam(bencher: &mut Bencher) {
        let lag: LaguerrePolynomial<f64> = LaguerrePolynomial::new(1_000);
        bencher.iter(|| lag.zeros_amsterdam())
    }
}
