mod problems;

use std::iter::Sum;

use integrator::adaptive_quadrature::simpson::adaptive_simpson_method;
use num::Float;
use problems::Problem;

pub fn adaptive_simpson_problems<F: Float + Send + Sum + Sync>() -> Vec<Problem<F>> {
    vec![
        problems::problem01(),
        problems::problem02(),
        problems::problem03(),
        problems::problem04(),
        problems::problem05(),
        problems::problem06(),
        problems::problem07(),
        problems::problem08(),
        problems::problem09(),
        problems::problem10(),
        problems::problem11(),
        problems::problem12(),
        problems::problem13(),
        problems::problem14(),
        problems::problem15(),
        problems::problem16(),
        problems::problem17(),
        problems::problem18(),
        problems::problem19(),
        problems::problem20(),
        problems::problem21(),
        problems::problem22(),
        problems::problem23(),
        problems::problem24(),
        problems::problem25(),
        problems::problem26(),
        problems::problem29(),
        problems::problem30(),
    ]
}

fn test_problem_f64(problem: Problem<f64>) {
    let f = problem.function;
    let (a, b) = problem.limits;

    let tolerance = 10.0e-6;
    let min_h = 10.0e-3;

    let result = adaptive_simpson_method(f, a, b, min_h, tolerance);

    match result {
        Ok(res) => {
            let test_passed = problem.check_result(res);
            let test_result = if test_passed { "passed" } else { "failed" };

            println!("{} -- {}", res, problem.exact);

            println!(
                "Method:AdaptiveSimpson -- Problem number:{} -- test:{}",
                problem.id, test_result
            );
            assert!(problem.check_result(res));
        }
        Err(err) => println!(
            "Method:AdaptiveSimpson -- Problem number:{} -- {}",
            problem.id, err
        ),
    };
}

fn test_problem_f32(problem: Problem<f32>) {
    let f = problem.function;
    let (a, b) = problem.limits;

    let tolerance = 10.0e-6;
    let min_h = 10.0e-3;

    let result = adaptive_simpson_method(f, a, b, min_h, tolerance);

    match result {
        Ok(res) => {
            let test_passed = problem.check_result(res);
            let test_result = if test_passed { "passed" } else { "failed" };

            println!("{} -- {}", res, problem.exact);

            println!(
                "Method:AdaptiveSimpson -- Problem number:{} -- test:{}",
                problem.id, test_result
            );
            assert!(problem.check_result(res));
        }
        Err(err) => println!(
            "Method:AdaptiveSimpson -- Problem number:{} -- {}",
            problem.id, err
        ),
    };
}

#[test]
fn test_f32_problems() {
    let problems: Vec<Problem<f32>> = adaptive_simpson_problems();

    for problem in problems.into_iter() {
        test_problem_f32(problem);
    }
}

#[test]
fn test_f64_problems() {
    let problems: Vec<Problem<f64>> = adaptive_simpson_problems();

    for problem in problems.into_iter() {
        test_problem_f64(problem);
    }
}
