use integrator::newton_cotes::{
    rectangle::rectangle_rule, simpson::simpson_rule, trapezoidal::trapezoidal_rule,
};

const A: f64 = 0.0;
const B: f64 = 1.0;
const NUM_STEPS: usize = 1_000_000_000;

fn rectangle() {
    fn square(x: f64) -> f64 {
        x.powi(2)
    }

    let integral = rectangle_rule(square, A, B, NUM_STEPS);
    println!("{}", integral)
}

fn newton() {
    fn square(x: f64) -> f64 {
        x.powi(2)
    }

    let integral = rectangle_rule(square, A, B, NUM_STEPS);
    println!("{}", integral)
}

fn trapezoidal() {
    fn square(x: f64) -> f64 {
        x.powi(2)
    }

    let integral = trapezoidal_rule(square, A, B, NUM_STEPS);
    println!("{}", integral)
}

fn simpson() {
    fn square(x: f64) -> f64 {
        x.powi(2)
    }

    let integral = simpson_rule(square, A, B, NUM_STEPS);
    println!("{}", integral)
}

fn main() {
    rectangle();
    trapezoidal();
    simpson();
    newton()
}
