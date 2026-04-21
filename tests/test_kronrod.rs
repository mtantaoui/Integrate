mod problems;

use integrate::gauss_kronrod::gauss_kronrod_rule;
use problems::{
    problem01, problem04, problem05, problem08, problem10, problem11, problem18, problem20,
    problem22, problem26, Problem,
};

// Order-7 Kronrod rule: 15-point evaluation on the whole interval.
const N: usize = 7;

/// Smooth, finite-limit problems suitable for a single-interval 15-point GK rule.
pub fn kronrod_problems_f64() -> Vec<Problem<f64>> {
    vec![
        problem01(), // ∫₀¹ eˣ dx
        problem04(), // ∫₋₁¹ 0.92 cosh x − cos x dx
        problem05(), // ∫₋₁¹ 1/(x⁴+x²+0.9) dx
        problem08(), // ∫₀¹ 1/(1+x⁴) dx
        problem10(), // ∫₀¹ 1/(1+x) dx
        problem11(), // ∫₀¹ 1/(1+eˣ) dx
        problem18(), // ∫₀¹ x/(eˣ+1) dx
        problem20(), // ∫₋₁¹ 1/(x²+1.005) dx
        problem22(), // ∫₀¹ 1/(x⁴+x²+1) dx
        problem26(), // ∫₀²π exp(cos x) dx
    ]
}

fn test_problem_f64(problem: Problem<f64>) {
    let f = problem.function;
    let (a, b) = problem.limits;

    let (integral, _error) = gauss_kronrod_rule(f, a, b, N).unwrap();

    let test_passed = problem.check_result(integral);
    let test_result = if test_passed { "passed" } else { "failed" };

    println!(
        "Method:GaussKronrod -- Problem number:{} -- test:{}",
        problem.id, test_result
    );

    assert!(
        test_passed,
        "Problem {}: got {integral}, expected {}",
        problem.id, problem.exact
    );
}

#[test]
fn test_f64_problems() {
    for problem in kronrod_problems_f64() {
        test_problem_f64(problem);
    }
}
