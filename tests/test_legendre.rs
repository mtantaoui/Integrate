mod problems;

const N: usize = 100_000;

use std::iter::Sum;

use integrate::gauss_quadrature::legendre::legendre_rule;
use num::Float;
use problems::{
    problem01, problem02, problem03, problem04, problem05, problem06, problem07, problem08,
    problem09, problem10, problem11, problem12, problem14, problem15, problem16, problem17,
    problem18, problem19, problem20, problem21, problem22, problem23, problem24, problem25,
    problem26, problem27, problem28, problem29, problem30, Problem,
};

pub fn legendre_problems<F: Float + Send + Sum + Sync>() -> Vec<Problem<F>> {
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
        problem27(),
        problem28(),
        problem29(),
        problem30(),
    ]
}

fn test_problem_f64(problem: Problem<f64>) {
    let f = problem.function;
    let (a, b) = problem.limits;
    let n: usize = N;

    let result = legendre_rule(f, a, b, n);

    let test_passed = problem.check_result(result);
    let test_result = if test_passed { "passed" } else { "failed" };

    println!(
        "Method:Legendre -- Problem number:{} -- test:{}",
        problem.id, test_result
    );

    assert!(problem.check_result(result));
}

fn test_problem_f32(problem: Problem<f32>) {
    let f = problem.function;
    let (a, b) = problem.limits;
    let n: usize = N;

    let result = legendre_rule(f, a, b, n) as f32;

    let test_passed = problem.check_result(result);
    let test_result = if test_passed { "passed" } else { "failed" };

    println!(
        "Method:Legendre -- Problem number:{} -- test:{}",
        problem.id, test_result
    );

    assert!(problem.check_result(result));
}

#[test]
fn test_f32_problems() {
    let problems: Vec<Problem<f32>> = legendre_problems();

    for problem in problems.into_iter() {
        test_problem_f32(problem);
    }
}

#[test]
fn test_f64_problems() {
    let problems: Vec<Problem<f64>> = legendre_problems();

    for problem in problems.into_iter() {
        test_problem_f64(problem);
    }
}
