use integrator::newton_cotes::rectangle::rectangle_rule;

fn main() {
    fn square(x: f64) -> f64 {
        x.sin()
    }

    let a = -1.0;
    let b = 1.0;

    let num_steps: usize = 1_000_000_000;

    let integral = rectangle_rule(square, a, b, num_steps);
    println!("{}", integral)
}
