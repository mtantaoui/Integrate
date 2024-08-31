use integrator::{
    gauss_quadrature::laguerre::gauss_laguerre_rule,
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

fn laguerre() {
    fn f(x: f64) -> f64 {
        (-x).exp()
    }

    let n = NUM_STEPS;
    let value = gauss_laguerre_rule(f, n);

    println!("laguerre: {}", value);
}

fn main() {
    laguerre();

    romberg();
    rectangle();
    trapezoidal();
    simpson();
    newton();
}
