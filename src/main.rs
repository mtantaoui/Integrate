use integrator::newton_cotes::{
    rectangle::rectangle_rule, simpson::simpson_rule, trapezoidal::trapezoidal_rule,
};

fn rectangle() {
    fn square(x: f64) -> f64 {
        x.powi(2)
    }

    let a = 0.0;
    let b = 1.0;

    let num_steps: usize = 1_000_000;

    let integral = rectangle_rule(square, a, b, num_steps);
    println!("{}", integral)
}

fn trapezoidal() {
    fn square(x: f64) -> f64 {
        x.powi(2)
    }

    let a = 0.0;
    let b = 1.0;

    let num_steps: usize = 1_000_000;

    let integral = trapezoidal_rule(square, a, b, num_steps);
    println!("{}", integral)
}

fn simpson() {
    fn square(x: f64) -> f64 {
        x.powi(2)
    }

    let a = 0.0;
    let b = 1.0;

    let num_steps: usize = 1_000_000_000;

    let integral = simpson_rule(square, a, b, num_steps);
    println!("{}", integral)
}

fn main() {
    rectangle();
    trapezoidal();
    simpson();
}
