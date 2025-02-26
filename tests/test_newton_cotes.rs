mod problems;

use std::iter::Sum;

use integrate::newton_cotes::{newton_rule, rectangle_rule, simpson_rule, trapezoidal_rule};
use num::Float;

use itertools::Itertools;
use problems::{
    problem01, problem02, problem03, problem04, problem05, problem06, problem07, problem08,
    problem09, problem10, problem11, problem12, problem13, problem14, problem15, problem16,
    problem17, problem18, problem19, problem20, problem21, problem22, problem23, problem24,
    problem25, problem26, problem27, problem28, problem29, problem30, Problem,
};

pub fn newton_cotes_problems<F: Float + Send + Sum + Sync>() -> Vec<Problem<F>> {
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
        problem27(),
        problem28(),
        problem29(),
        problem30(),
    ]
}

#[derive(Debug, Clone, Copy)]
enum Methods {
    Rectangle,
    Trapezoidal,
    Newton3Over8,
    Simpson,
}

impl Methods {
    pub fn iter() -> core::array::IntoIter<Methods, 4> {
        [
            Methods::Rectangle,
            Methods::Trapezoidal,
            Methods::Newton3Over8,
            Methods::Simpson,
        ]
        .into_iter()
    }

    pub fn display(&self) -> &str {
        match self {
            Methods::Rectangle => "Rectangle rule",
            Methods::Trapezoidal => "Trapezoidal Rule",
            Methods::Newton3Over8 => "Newton 3/8 Rule",
            Methods::Simpson => "Simpson Method",
        }
    }
}

fn test_problem_f64(problem: Problem<f64>, method: Methods) {
    let f = problem.function;
    let (a, b) = problem.limits;
    let n: usize = problem.n;

    let result = integrate(method, f, a, b, n);

    let test_passed = problem.check_result(result);
    let test_result = if test_passed { "passed" } else { "failed" };

    println!(
        "Method:{} -- Problem number:{} -- test:{}",
        method.display(),
        problem.id,
        test_result
    );

    assert!(problem.check_result(result));
}

fn test_problem_f32(problem: Problem<f32>, method: Methods) {
    let f = problem.function;
    let (a, b) = problem.limits;
    let n: usize = problem.n;

    let result = integrate(method, f, a, b, n) as f32;

    let test_passed = problem.check_result(result);
    let test_result = if test_passed { "passed" } else { "failed" };

    println!(
        "Method:{} -- Problem number:{} -- test:{}",
        method.display(),
        problem.id,
        test_result
    );

    assert!(problem.check_result(result));
}

fn integrate<F: Float + Send + Sync>(method: Methods, f: fn(F) -> F, a: F, b: F, n: usize) -> f64 {
    match method {
        Methods::Rectangle => rectangle_rule(f, a, b, n),
        Methods::Trapezoidal => trapezoidal_rule(f, a, b, n),
        Methods::Newton3Over8 => newton_rule(f, a, b, n),
        Methods::Simpson => simpson_rule(f, a, b, n),
    }
}

#[test]
fn test_f32_problems() {
    let problems: Vec<Problem<f32>> = newton_cotes_problems();
    let methods = Methods::iter();

    for (problem, method) in problems.into_iter().cartesian_product(methods.into_iter()) {
        test_problem_f32(problem, method);
    }
}

#[test]
fn test_f64_problems() {
    let problems: Vec<Problem<f64>> = newton_cotes_problems();
    let methods = Methods::iter();

    for (problem, method) in problems.into_iter().cartesian_product(methods.into_iter()) {
        test_problem_f64(problem, method);
    }
}
