use std::f64::consts::PI;

use integrator::{
    adaptive_quadrature::simpson::adaptive_simpson_method,
    gauss_quadrature::{
        chebyshev::{gauss_first_kind_chebyshev_rule, gauss_second_kind_chebyshev_rule},
        hermite::gauss_hermite_rule,
        laguerre::gauss_laguerre_rule,
    },
    newton_cotes::{
        rectangle::rectangle_rule, simpson::simpson_rule, trapezoidal::trapezoidal_rule,
    },
    romberg::romberg_method,
};

const A: f64 = 0.0;
const B: f64 = 1.0;
const NUM_STEPS: usize = 300;
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
    let n: usize = 20;
    let integral = romberg_method(square, A, B, n);
    println!("romberg: {}", integral)
}

fn hermite() {
    fn f(_: f64) -> f64 {
        1.0
    }

    let integral = gauss_laguerre_rule(f, 170); // 170  is the maximum value for n, to test on x86 arch
    println!("hermite integral {}\t value {}", integral, 1);

    let integral = gauss_hermite_rule(f, 170); // 170  is the maximum value for n, to test on x86 arch
    println!("hermite integral {}\t value {}", integral, PI.sqrt());
}

fn chebyshev() {
    fn f(x: f64) -> f64 {
        (1.0 - x.powi(2)).sqrt()
    }

    fn g(x: f64) -> f64 {
        1.0 / (1.0 - x.powi(2)).sqrt()
    }

    let n: usize = 100;

    let integral = gauss_first_kind_chebyshev_rule(f, n);
    println!("integral {}\t value {}", integral, 2);

    let integral = gauss_second_kind_chebyshev_rule(g, n);
    println!("integral {}\t value {}", integral, 2);
}

fn adaptive_simpson() {
    fn square(x: f64) -> f64 {
        x.powi(POW)
    }

    let tolerance = 10.0e-6;
    let min_h = 10.0e-3;

    let integral = adaptive_simpson_method(square, A, B, min_h, tolerance);
    println!("adaptive simpson: {}", integral.unwrap());
}

fn main() {
    adaptive_simpson();
    chebyshev();
    hermite();
    romberg();
    rectangle();
    trapezoidal();
    simpson();
    newton();
}
