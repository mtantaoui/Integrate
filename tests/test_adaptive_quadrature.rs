mod problems;

use std::iter::Sum;

use integrate::adaptive_quadrature::simpson::adaptive_simpson_method;
use num::Float;

use problems::{
    problem01, problem02, problem03, problem04, problem05, problem06, problem07, problem08,
    problem09, problem10, problem11, problem12, problem13, problem14, problem15, problem16,
    problem17, problem18, problem19, problem20, problem21, problem22, problem23, problem24,
    problem25, problem26, problem29, problem30, Problem,
};

pub fn adaptive_simpson_problems<F: Float + Send + Sum + Sync>() -> Vec<Problem<F>> {
    vec![
        problem01(),
        problem02(),
        problem03(),
        problem04(),
        problem05(),
        problem06(),
        problem07(),
        problem08(),
        problem09(),
        problem10(),
        problem11(),
        problem12(),
        problem13(),
        problem14(),
        problem15(),
        problem16(),
        problem17(),
        problem18(),
        problem19(),
        problem20(),
        problem21(),
        problem22(),
        problem23(),
        problem24(),
        problem25(),
        problem26(),
        problem29(),
        problem30(),
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
