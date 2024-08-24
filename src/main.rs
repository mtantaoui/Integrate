use integrator::{
    newton_cotes::{
        rectangle::rectangle_rule, simpson::simpson_rule, trapezoidal::trapezoidal_rule,
    },
    romberg::romberg_method,
    utils::orthogonal_polynomials::{LaguerrePolynomial, OrthogonalPolynomial},
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
    let lag: LaguerrePolynomial<f64> = LaguerrePolynomial::new(n);
    let roots = lag.zeros();
    let ams_roots = lag.zeros_amsterdam();
    let ub = lag.zeros_upper_bounds();
    let lb = lag.zeros_lower_bounds();

    for i in 0..n {
        println!("{} ", lag.eval(roots[i]))
    }

    // for (i, ((root, upper_bound), lower_bound)) in roots.into_iter().zip(ub).zip(lb).enumerate() {
    //     println!(
    //         "{} - root:{} - upper:{} - lower:{} - {}",
    //         i,
    //         root,
    //         upper_bound,
    //         lower_bound,
    //         root <= upper_bound
    //     )
    // }
}

fn main() {
    laguerre_roots();

    romberg();
    rectangle();
    trapezoidal();
    simpson();
    newton();
}
