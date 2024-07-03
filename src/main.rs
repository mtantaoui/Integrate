use integrator::{
    gauss_quadrature::laguerre::eigenvalue,
    newton_cotes::{
        rectangle::rectangle_rule, simpson::simpson_rule, trapezoidal::trapezoidal_rule,
    },
    romberg::romberg_method,
};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

const A: f64 = 0.0;
const B: f64 = 1.0;
const NUM_STEPS: usize = 20;
const POW: i32 = 2;

fn rectangle() {
    fn square(x: f64) -> f64 {
        x.powi(POW)
    }

    let integral = rectangle_rule(square, A, B, NUM_STEPS);
    println!("rectangle: {}", integral)
}

fn newton() {
    fn square(x: f64) -> f64 {
        x.powi(POW)
    }

    let integral = rectangle_rule(square, A, B, NUM_STEPS);
    println!("newton: {}", integral)
}

fn trapezoidal() {
    fn square(x: f64) -> f64 {
        x.powi(POW)
    }

    let integral = trapezoidal_rule(square, A, B, NUM_STEPS);
    println!("trapezoidal: {}", integral)
}

fn simpson() {
    fn square(x: f64) -> f64 {
        x.powi(POW)
    }

    let integral = simpson_rule(square, A, B, NUM_STEPS);
    println!("simpson: {}", integral)
}

fn romberg() {
    fn square(x: f64) -> f64 {
        x.powi(POW)
    }

    let integral = romberg_method(square, A, B, NUM_STEPS);
    println!("romberg: {}", integral)
}

fn givens_test() {
    // let diagonal = vec![1.0, 2.0, 3.0, 4.0];
    // let mut offdiagonal = vec![0.0, 1.0, 2.0, 3.0];

    let diagonal = vec![2.0, 1.0, 0.0, 4.0];
    let offdiagonal = vec![0.0, 3.0, 0.0, 1.0];

    // let diagonal = vec![1.0, 1.0, 1.0, 1.0];
    // let mut offdiagonal = vec![0.0, 0.0, 0.0, 0.0];

    let n = diagonal.len();

    let eigenvalues: Vec<f64> = (0..n)
        .into_par_iter()
        .map(|k| eigenvalue(diagonal.as_slice(), offdiagonal.as_slice(), k))
        .collect();

    println!("eig: {:?}", eigenvalues);
}

fn main() {
    givens_test();
    romberg();
    rectangle();
    trapezoidal();
    simpson();
    newton();
}
