# Romberg's method

## Example

```rust,editable
use integrate::romberg::romberg_method;


fn square(x: f64) -> f64 {
    x.powi(2)
}

let a = 0.0;
let b = 1.0;

let num_steps: usize = 10;

let integral = romberg_method(square, a, b, num_steps);
println!("{}",integral);
```

## Understanding Romberg's method

Romberg's method is used to estimate the integral of a function on a closed and bounded interval.

Classically, the method consists of successively applying the composite trapezoidal rule, each time
halving the length of the subintervals, and using a linear combination of the resulting sequence of
estimates to estimate the integral, by successively deleting the low order error terms in
the Euler-Maclaurin summation formula.

The process terminates when the change of the estimate is within a preassigned tolerance, within
a preassigned number of successive estimates.

The Euler-Maclaurin summation formula relates the integral of a function \\(f(x)\\) over
a closed and bounded interval \\(\[a,b\]\\) , \\(\int\_{a}^{b} f(x) dx\\), and the composite trapezoidal rule,

\\[
T_h (f) = h \left[ \frac{f(a)}{2} + f(a+h) + ··· + f(b-h) + \frac{f(b)}{2} \right]
\\]

by

\\[
\begin{split}
T_h(f) &= \int\_{a}^{b} f(x) dx + \left(\frac{h^2}{12}\right) \left[f^\prime(b) - f^\prime(a)\right] - \left(\frac{h^4}{720}\right) \left[f^{(3)}(b) - f^{(3)}(a) \right] \\\\
&+ ... + K h^{2p-2} \left[ f^{(2p-3)}(b) - f^{(2p-3)(a)} \right] + O(h^{2p})
\end{split}
\\]

where \\(f^\prime\\), \\(f^{(3)}\\), and \\(f^{(2p-3)}\\) are the first, third and \\((p-3)rd\\) derivatives
of \\(f\\) and \\(K\\) is a constant.

If the subinterval length is halved, then

\\[
\begin{split}
T\_{\frac{h}{2}}(f) &= \int\_{a}^{b} f(x) dx + \frac{h^2}{4·12} \left[ f^\prime(b) - f^\prime(a) \right] - \left( \frac{h^4}{16·720} \right) \left[ f^{(3)}(b) - f^{(3)}(a) \right] \\\\
&+ ... + K \left( \frac{h}{2} \right)^{2p-2} \left[ f^{(2p-3)}(b) - f^{(2p-3)(a)} \right] + O(h^{2p})
\end{split}
\\]

So that

\\[
\begin{split}
\frac{4 T\_{\frac{h}{2}}(f) - T\_{h}(f)}{3} &= \int\_{a}^{b} f(x) dx + \left( \frac{h^4}{2880} \right) \left[ f^{(3)}(b) - f^{(3)}(a) \right] - \left( \frac{h^6}{96768} \right) \left[ f^{(5)}(b) - f^{(5)}(a) \right] \\\\
&+ ... + K h^{2p-2} \left[ f^{(2p-3)}(b) - f^{(2p-3)(a)} \right] + O(h^{2p})
\end{split}
\\]

The \\(h^2\\) term has vanished. This process can be continued, each halving of the subinterval
length results in a new composite trapezoidal rule estimate of the integral, which can be
combined with previous estimates to yield an estimate, in which the lowest order term
involving \\(h\\) vanishes.

The easiest way to combine all the estimates from applications
of the trapezoidal rule by halving the length of the subintervals is to arrange the estimates
in a column, from the coarsest estimate to the finest estimate. To form the second,
column take two adjacent values from the first column, subtract the finer estimate
from the coarser estimate divide by 3 and add to the finer estimate.

The second column is automatically arranged from the coarsest estimate to the finest
with one less element. Continue, form the third column by taking two adjacent values
from the second column, subtract the finer estimate from the coarser estimate,
divide by 15 and add to the finer estimate.

The third column is automatically arranged from the coarsest to the finest estimate with one fewer element than in the
second column.

This process is continued until there is only one element in the last column, this
is the estimate of the integral.

The numbers which are used the divide the difference of two adjacent elements in the \\(i^{th}\\) column is \\(4^i - 1\\).
