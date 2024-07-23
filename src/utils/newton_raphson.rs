use num::Float;

pub fn newton_raphson<F: Float, Function>(f: Function, df: Function, a: F, tolerance: F) -> F
where
    Function: Fn(F) -> F,
{
    // let mut a = a;
    let mut delta;
    loop {
        let fa = f(a);

        if fa.is_zero() {
            return a;
        }

        let dfa = df(a);

        // If an attempt to divide by zero, TODO: return a specific error
        if dfa.is_zero() {
            panic!()
        }

        delta = -fa / dfa;

        // If the location of the root relative to a is less than //
        // the tolerance, return the estimate of the root.        //

        if delta.abs() < tolerance {
            break;
        }
    }

    return a + delta;
}
