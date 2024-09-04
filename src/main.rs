use std::f64::consts::PI;

use integrator::{
    gauss_quadrature::{hermite::gauss_hermite_rule, laguerre::gauss_laguerre_rule},
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
    println!("integral {}\t value {}", integral, 1);

    let integral = gauss_hermite_rule(f, 170); // 170  is the maximum value for n, to test on x86 arch
    println!("integral {}\t value {}", integral, PI.sqrt());
}

fn main() {
    hermite();

    romberg();
    rectangle();
    trapezoidal();
    simpson();
    newton();
}
