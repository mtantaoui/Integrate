//! The rectangle rule approximates the integral of a function $f(x)$ on the
//! closed and bounded interval $[a, a+h]$ of length $h > 0$ by the (signed) area
//! of the rectangle with length h and height the value of the function $f(x)$
//! evaluated at the midpoint of the interval, $f(a+h/2)$.
//!
//! The composite rectangle rule is used to approximate the integral of a function
//! $f(x)$ over a closed and bounded interval $[a, b]$ where $a < b$, by decomposing
//! the interval $[a, b]$ into $n > 1$ subintervals of equal length $h = \frac{b-a}{n}$
//! and adding the results of applying the rectangle rule to each subinterval.
//!
//! By abuse of language both the composite rectangle rule and the rectangle rule sometimes
//! are referred to simply as the rectangle rule.
//!
//! Let $\int_{a}^{b} f(x) dx$ be the integral of f(x) over the closed and bounded interval $\[a ,b \]$,
//! and let $R_h(f)$ be the result of applying the rectangle rule with n subintervals of length h, i.e.
//! $$ R_h(f)=h [ f(a+h/2) + f(a+3h/2) + ··· + f(b-h/2) ] $$
//! An immediate consequence of the Euler-Maclaurin summation formula yields the following equation
//! relating $\int_{a}^{b} f(x) dx$ and $R_h(f)$
//! $$ R_h(f) = \int_{a}^{b} f(x) dx - \frac{h^2}{24} (f' (b) - f' (a) ) +  \frac{7h^4}{5760} ( f^{(3)}(b) - f^{(3)}(a) ) + $$
//! $$ ··· + K h^{2p-2} (f^{(2p-3)}(b) - f^{(2p-3)}(a) ) + O(h^{2p})  $$
//!
//! where $f'$, $f^{(3)}$, and $f^{(2p-3)}$ are the first, third and $(2p-3)rd$ derivatives of $f$ and $K$ is a constant.
//!
//! The last term, $O(h^{2p})$ is important. Given an infinitely differentiable function
//! in which the first $2p-3$ derivatives vanish at both endpoints of the interval of integration,
//! it is not true that $R_h(f) = \int_{a}^{b} f(x) dx$, but rather what the theorem says is that
//! $$ \lim_{h \to 0} \mid \frac{R_h(f) - \int_{a}^{b} f(x)dx}{h^{2p}} \mid < M $$
//! where $M>0$
//!
//! If $f$ is at least twice differentiable on the interval $\[a,b\]$, then applying the mean-value
//! theorem to
//! $$ R_h(f) - \int_{a}^{b} f(x) dx = -\frac{h^2}{24} (f' (b) - f' (a)) +  \frac{7h^4}{5760} \[ f^{(3)}(b) - f^{(3)}(a) \]+$$
//! $$ ··· + K h^{2p-2} (f^{(2p-3)}(b) - f^{(2p-3)}(a)) + O(h^{2p}) $$
//! yields the standard truncation error expression
//!
//! $$ R_h(f) - \int_{a}^{b} f(x) dx = -\frac{h^2}{24} (b - a) f''(c) $$
//! for some point $c$ where $a ≤ c ≤ b$.
//!
//! A corollary of which is that if $f''(x) = 0$ for all $x$ in $\[a,b\]$,
//! i.e. if $f(x)$ is linear, then the rectangle rule is exact.
//!
//! The Euler-Maclaurin summation formula also shows that usually $n$ should be chosen large enough
//! so that $h = (b - a) / n < 1$. For example, if h = 0.1 then
//! $$ R_{0.1}(f) = \int_{a}^{b} f(x) dx  - 0.00042 (f'(b) - f'(a)) + 0.00000012 (f''(b) - f''(a)) + ...   $$
//! and if $h = 0.01$ then
//! $$ R_{0.01}(f) = \int_{a}^{b} f(x) dx  - 0.0000042 (f'(b) - f'(a)) + 0.000000000012 (f''(b) - f''(a)) + ...   $$
//! while if $h=10$ then
//! $$ R_{0.01}(f) = \int_{a}^{b} f(x) dx  - 4.1667 (f'(b) - f'(a)) + 12.15 (f''(b) - f''(a)) + ... $$
//! However, if the function $f(x)$ is linear, then $n$ may be chosen to be $1$.

use num::{ToPrimitive, Unsigned};
use num_traits::real::Real;
use rayon::prelude::*;

/// This function integrates $f(x)$ from $a$ to $a+n\star h$ using the rectangle
///  rule by summing from the left end of the interval to the right end.
/// # Examples
///
/// ```
/// println!("This is Rectangle rule example placeholder");
///
///
/// ```
///
/// # Resources
/// [Methods of numerical Integration (2nd edition), by Philip J. Davis and Philip Rabinowitz.](https://www.cambridge.org/core/journals/mathematical-gazette/article/abs/methods-of-numerical-integration-2nd-edition-by-philip-j-davis-and-philip-rabinowitz-pp-612-3650-1984-isbn-0122063600-academic-press/C331158D0392E1D5CD9B0C6ED4EE5F43)
pub fn rectangle_rule<
    R1: Real + Sync,
    R2: Real + std::marker::Send,
    U: Unsigned + ToPrimitive + Copy,
>(
    f: fn(R1) -> R2,
    a: R1,
    b: R1,
    n: U,
) -> f64
// TODO: complete function documentation with example and more details
{
    let h: R1 = (b - a) / R1::from(n).unwrap();

    let integral: f64 = (0..(n.to_usize().unwrap()))
        .into_par_iter()
        .map(|i| {
            let i = R1::from(i).unwrap();
            let x = a + i * h + (h / R1::from(2).unwrap());
            f(x).to_f64().unwrap()
        })
        .sum();
    integral * h.to_f64().unwrap()
}

#[cfg(test)]
mod tests {
    use std::ops::Div;

    use super::*;
    // TODO: find a way to compare two approxiamte floats in rust (almost_equal in numpy equivalent)
    #[test]
    fn test_integral_value() -> Result<(), String> {
        fn f1(x: f64) -> f64 {
            x.powi(2)
        }

        let a = 0.0;
        let b = 1.0;

        let num_steps: usize = 1_000_000;

        let integral = rectangle_rule(f1, a, b, num_steps);
        let epsilon = 10e-10;

        assert!((integral - 1.0.div(3.0) as f64).abs() < epsilon);

        Ok(())
    }
}
