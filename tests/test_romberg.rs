#[path = "./problems.rs"]
mod pbs;
const N: usize = 20;

use pbs::{problems_vec, Problem};

use integrator::romberg::romberg_method;

fn test_problem_f64(problem: Problem<f64>) {
    let f = problem.function;
    let (a, b) = problem.limits;
    let n: usize = N;

    let result = romberg_method(f, a, b, n);

    let test_passed = problem.check_result(result);
    let test_result = if test_passed { "passed" } else { "failed" };

    println!(
        "Method:Romberg -- Problem number:{} -- test:{}",
        problem.id, test_result
    );

    assert!(problem.check_result(result));
}

fn test_problem_f32(problem: Problem<f32>) {
    let f = problem.function;
    let (a, b) = problem.limits;
    let n: usize = N;

    let result = romberg_method(f, a, b, n) as f32;

    let test_passed = problem.check_result(result);
    let test_result = if test_passed { "passed" } else { "failed" };

    println!(
        "Method:Romberg -- Problem number:{} -- test:{}",
        problem.id, test_result
    );

    assert!(problem.check_result(result));
}

#[test]
fn test_f32_problems() {
    let problems: Vec<Problem<f32>> = problems_vec();

    for problem in problems.into_iter() {
        test_problem_f32(problem);
    }
}

#[test]
fn test_f64_problems() {
    let problems: Vec<Problem<f64>> = problems_vec();

    for problem in problems.into_iter() {
        test_problem_f64(problem);
    }
}
