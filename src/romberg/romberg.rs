//! The Euler-Maclaurin summation formula relates the integral of a function $f(x)$ over
//! a closed and bounded interval $\[a,b\]$ , $\int_{a}^{b} f(x) dx$, and the composite trapezoidal rule,
//!
//! $$ T_h(f) = h \[ f(a)/2 + f(a+h) + ··· + f(b-h) + f(b)/2 \]$$
//!
//! by
//!
//! $$ T_h(f) = \int_{a}^{b} f(x) dx + \frac{h^2}{12}\[f'(b) - f'(a)\] - (\frac{h^4}{720})[f^{(3)}(b) - f^{(3)}(a)]$$
//! $$ + ... + K h^{2p-2} \[f^{(2p-3)}(b) - f^{(2p-3)(a)}\] + O(h^{2p})$$
//!
//! where $f'$, $f^{(3)}$, and $f^{(2p-3)}$ are the first, third and $(p-3)rd$ derivatives
//! of $f$ and $K$ is a constant.
//!
//! If the subinterval length is halved, then
//!
//! $$ T_{\frac{h}{2}}(f) = \int_{a}^{b} f(x) dx + \frac{h^2}{4·12}\[f'(b) - f'(a)\] - (\frac{h^4}{16·720})[f^{(3)}(b) - f^{(3)}(a)]$$
//! $$ + ... + K (\frac{h}{2})^{2p-2} \[f^{(2p-3)}(b) - f^{(2p-3)(a)}\] + O(h^{2p})$$
//!
//! So that
//!
//! $$ \frac{4T_{\frac{h}{2}}(f) - T_{h}(f)}{3} = \int_{a}^{b} f(x) dx + \frac{h^4}{2880}\[f^{(3)}(b) - f^{(3)}(a)\] - (\frac{h^6}{96768})[f^{(5)}(b) - f^{(5)}(a)]$$
//! $$ + ... + K h^{2p-2} \[f^{(2p-3)}(b) - f^{(2p-3)(a)}\] + O(h^{2p})$$
//!
//! The $h^2$ term has vanished. This process can be continued, each halving of the subinterval
//! length results in a new composite trapezoidal rule estimate of the integral which can be
//! combined with previous estimates to yield an estimate in which the lowest order term
//! involving $h$ vanishes. The easiest way to combine all the estimates from applications
//! of the trapezoidal rule by halving the length of the subintervals is to arrange the estimates
//! in a column from the coarsest estimate to the finest estimate. To form the second,
//! column take two adjacent values from the first column, subtract the finer estimate
//! from the coarser estimate divide by 3 and add to the finer estimate.
//!
//! The second column is automatically arranged from the coarsest estimate to the finest
//! with one less element. Continue, form the third column by taking two adjacent values
//! from the second column, subtract the finer estimate from the coarser estimate,
//! divide by 15 and add to the finer estimate. The third column is automatically
//! arranged from the coarsest to the finest estimate with one fewer element than in the
//! second column.
//!
//! This process is continued until there is only one element in the last column, this
//! is the estimate of the integral.
//!
//! The numbers which are used the divide the difference of two adjacent elements in the $i^{th}$ column is $4^i - 1$.

use std::cmp::max;

use num::{Float, ToPrimitive, Unsigned};

/// If $T_h(f)$ is the result of applying the trapezoidal rule to approximating
/// the integral of $f(x)$ on $\[a, b\]$ using subintervals of length $h$,   
/// then if $\int_{a}^{b} f(x) dx$ is the integral of $f(x)$ on $\[a,b\]$, then                 
///                           \int_{a}^{b} f(x) dx = \lim_{h \to 0}  T_h(f)                       
/// where the limit is taken as h approaches 0.    
///                          
/// The classical Romberg method applies Richardson Extrapolation to the    
/// limit of the sequence T_h(f), T_{\frac{h}{2}}(f), T_{\frac{h}{2}}(f), ... ,     
/// in which the limit is approached by successively deleting error terms   
/// in the Euler-MacLaurin summation formula.  
///
/// # Examples
/// ```
/// use integrator::romberg::romberg_method;
///
///
/// fn square(x: f64) -> f64 {
///     x.powi(2)
/// }
///
/// let a = 0.0;
/// let b = 1.0;
///
/// let num_steps: usize = 1_000_000;
///
/// let integral = newton_rule(square, a, b, num_steps);
/// ```
///
/// # Resources
/// [Methods of numerical Integration (2nd edition), by Philip J. Davis and Philip Rabinowitz.](https://www.cambridge.org/core/journals/mathematical-gazette/article/abs/methods-of-numerical-integration-2nd-edition-by-philip-j-davis-and-philip-rabinowitz-pp-612-3650-1984-isbn-0122063600-academic-press/C331158D0392E1D5CD9B0C6ED4EE5F43)
pub fn romberg_method<
    F1: Float + Sync + Ord,
    F2: Float + Send,
    U: Unsigned + ToPrimitive + Copy,
>(
    f: fn(F1) -> F2,
    a: F1,
    b: F1,
    n: U,
) -> f64 {
    // first term of trapezoid rule estimation
    let i_0 = (f(a) + f(b)).div(F2::from(2).unwrap()).to_f64().unwrap();

    todo!();

    // Initialize err and the first trapezoid estimates column, dt[0], to the numerical estimate
    // of the integral using the trapezoidal rule with a step size of h.

    Float::max(a, b).to_f64().unwrap()
}
