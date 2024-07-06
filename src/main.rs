use integrator::{
    gauss_quadrature::laguerre::{
        laguerre_polynomial_approximate_zeros, laguerre_polynomial_zeros,
    },
    newton_cotes::{
        rectangle::rectangle_rule, simpson::simpson_rule, trapezoidal::trapezoidal_rule,
    },
    romberg::romberg_method,
};

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

fn laguerre_roots() {
    let n = 100;
    let computed = laguerre_polynomial_zeros(n);
    let mut approximated = laguerre_polynomial_approximate_zeros(n);

    approximated.reverse();

    for (c, a) in computed.into_iter().zip(approximated) {
        println!("computed: {}", c);
        println!("approximated: {}", a);

        println!();
        println!();
    }
}

fn main() {
    laguerre_roots();

    romberg();
    rectangle();
    trapezoidal();
    simpson();
    newton();
}
