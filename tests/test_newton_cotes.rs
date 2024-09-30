#[path = "./problems.rs"]
mod pbs;

use num::Float;
use pbs::{problems_vec, Methods, Problem};

use integrator::newton_cotes::{
    newton::newton_rule, rectangle::rectangle_rule, simpson::simpson_rule,
    trapezoidal::trapezoidal_rule,
};

use itertools::Itertools;

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
    let problems: Vec<Problem<f32>> = problems_vec();
    let methods = Methods::iter();

    for (problem, method) in problems.into_iter().cartesian_product(methods.into_iter()) {
        test_problem_f32(problem, method);
    }
}

#[test]
fn test_f64_problems() {
    let problems: Vec<Problem<f64>> = problems_vec();
    let methods = Methods::iter();

    for (problem, method) in problems.into_iter().cartesian_product(methods.into_iter()) {
        test_problem_f64(problem, method);
    }
}
