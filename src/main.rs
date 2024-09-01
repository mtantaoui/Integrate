use integrator::{
    gauss_quadrature::hermite::Hermite,
    newton_cotes::{
        rectangle::rectangle_rule, simpson::simpson_rule, trapezoidal::trapezoidal_rule,
    },
    romberg::romberg_method,
    utils::orthogonal_polynomials::OrthogonalPolynomial,
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
    let h1: Hermite<f64> = Hermite::new(1);
    let h1_zeros = h1.zeros();

    let h2: Hermite<f64> = Hermite::new(2);
    let h2_zeros = h2.zeros();

    let h3: Hermite<f64> = Hermite::new(3);
    let h3_zeros = h3.zeros();

    let h4: Hermite<f32> = Hermite::new(4);
    let h4_zeros = h4.zeros();

    let h5: Hermite<f32> = Hermite::new(5);
    let h5_zeros = h5.zeros();

    println!("{:?}  \n", h1_zeros);
    println!("{:?}  \n", h2_zeros);
    println!("{:?}  \n", h3_zeros);
    println!("{:?}  \n", h4_zeros);
    println!("{:?}  \n", h5_zeros);
}

fn main() {
    laguerre();

    romberg();
    rectangle();
    trapezoidal();
    simpson();
    newton();
}
