# Newton's 3/8 rule

## Motivation

Newton's 3/8 rule fits a cubic polynomial through four equally-spaced points per subinterval and integrates it exactly. It achieves the same order of accuracy as Simpson's rule (\\(O(h^4)\\)) but uses three subintervals at a time, making it useful when the number of subintervals is a multiple of three.

## Example

```rust, editable
use integrate::newton_cotes::newton_rule;


let square = |x: f64| x * x;

let a = 0.0;
let b = 1.0;

let num_steps: usize = 1_000_000;

let integral = newton_rule(square, a, b, num_steps);
```

## Understanding Newton's 3/8 rule

Newton's \\(3/8\\) rule approximates the integral of a function \\(f(x)\\) on the closed and bounded
interval \\(\[a, a+h\]\\) of length \\(h > 0\\) by the integral on \\(\[a, a+h\]\\) of the cubic passing
through the points \\((a, f(a))\\), \\((a+\frac{h}{3}, f(a+\frac{h}{3}))\\), \\((a+\frac{2h}{3}, f(a+\frac{2h}{3}))\\) and \\((a+h, f(a+h))\\).

The composite Newton's 3/8 rule is used to approximate the integral of a function \\(f(x)\\)
over a closed and bounded interval \\(\[a, b\]\\) where \\(a < b\\), by decomposing the interval
\\(\[a, b\]\\) into \\(n > 1\\) subintervals of equal length \\(h = \dfrac{b-a}{n}\\), then adding the
results of applying the Newton's 3/8 rule to each subinterval.

By abuse of language both the composite Newton's 3/8 rule and Newton's 3/8
rule are referred to simply as Newton's 3/8 rule. Let \\(\int\_{a}^{b} f(x)dx\\) be the
integral of \\(f(x)\\) over the closed and bounded interval \\(\[a, b\]\\), and let \\(N_h(f)\\)
be the result of applying the Newton's 3/8 rule with \\(n\\) subintervals of length \\(h\\), i.e.

\\[
\begin{align}
N_h(f) &= \frac{h}{8} \left[ f(a) + 3 f\left(a+ \frac{h}{3} \right) + 3f\left(a+ \frac{2h}{3} \right) + 2 f(a + h) \right.\\\\
& \left. + \cdots + 2f(b-h) + 3f \left( b-\frac{2h}{3} \right) + 3f \left( b - \frac{h}{3} \right) + f(b) \right]
\end{align}
\\]

An immediate consequence of the Euler-Maclaurin summation formula relates \\(\int\_{a}^{b} f(x)dx\\) and \\(N_h(f)\\)

\\[
\begin{align}
N_h(f) &= \int\_{a}^{b} f(x)dx + \frac{h^4}{6480} \left[ f^{(3)}(b) - f^{(3)}(a) \right] - \frac{h^6}{244944} \left[ f^{(5)}(b) - f^{(5)}(a) \right] \\\\
& + \cdots + K h^{2p-2} \left[ f^{(2p - 3)}(b) - f^{(2p - 3)}(a) \right] + O(h^{2p})
\end{align}
\\]

where \\(f^{(3)}\\), \\(f^{(5)}\\), and \\(f^{(2p-3)}\\) are the third, fifth and \\((p-3)rd\\) derivatives
of \\(f\\) and \\(K\\) is a constant.

The last term, \\(O(h^{2p})\\) is important. Given an infinitely differentiable function in which the
first \\(2p-3\\) derivatives vanish at both endpoints of the interval of integration, it is not true that

\\[
N_h(f) = \int\_{a}^{b} f(x) dx
\\]

but rather what the theorem says is that

\\[
\lim_{h \to 0} \left| \frac{N_h(f) - \int_{a}^{b} f(x)dx}{h^{2p}} \right| < M
\\]

where \\(M > 0\\).

If \\(f\\) is at least four times differentiable on the interval \\(\[a,b\]\\), then applying the
mean-value theorem to

\\[
\begin{align}
N_h(f) - \int\_{a}^{b} f(x)dx &= \frac{h^4}{6480} \left[ f^{(3)}(b) - f^{(3)}(a) \right] - \frac{h^6}{244944} \left[ f^{(5)}(b) - f^{(5)}(a) \right] \\\\
& + \cdots + K h^{2p - 2} \left[ f^{(2p - 3)}(b) - f^{(2p - 3)}(a) \right] + O(h^{2p})
\end{align}
\\]

yields the standard truncation error expression

\\[
N_h(f) - \int_{a}^{b} f(x)dx = \frac{h^4}{6480} (b-a) f^{(4)}(c)
\\]

for some point \\(c\\) where \\(a \leq c \leq b\\).

A corollary of which is that if \\(f^{(4)}(x) = 0\\) for all \\(x\\) in \\(\[a,b\]\\), i.e. if \\(f(x)\\) is a cubic,
then Newton's 3/8 rule is exact.

The Euler-Maclaurin summation formula also shows that usually \\(n\\) should be chosen large enough so that \\(h = \dfrac{b-a}{n}< 1\\).

For example, if \\(h = 0.1\\) then

\\[
\begin{align}
N\_{0.1}(f) &= \int\_{a}^{b} f(x)dx + 1.5·10^{-8} \left[ f^{(3)}(b) - f^{(3)}(a) \right] \\\\
&- 4.1 \cdot 10^{-12} \left[ f^{(5)}(b) - f^{(5)}(a) \right] + \cdots
\end{align}
\\]

and if \\(h = 0.01\\) then

\\[
\begin{align}
N\_{0.01}(f) &= \int\_{a}^{b} f(x)dx + 1.5 \cdot 10^{-12} \left[ f^{(3)}(b) - f^{(3)}(a) \right] \\\\
&- 4.1 \cdot 10^{-18} \left[ f^{(5)}(b) - f^{(5)}(a) \right] + \cdots
\end{align}
\\]

while if \\(h = 10\\) then

\\[
\begin{align}
N\_{10}(f) &= \int\_{a}^{b} f(x)dx + 1.54 \left[ f^{(3)}(b) - f^{(3)}(a) \right] \\\\
&- 4.08 \left[ f^{(5)}(b) - f^{(5)}(a) \right] + \cdots
\end{align}
\\]

However, if the function \\(f(x)\\) is a cubic, then \\(n\\) may be chosen to be \\(1\\).

## Limitations

Newton's 3/8 rule converges as \\(O(h^4)\\), the same order as Simpson's rule, but requires the number of subintervals to be a multiple of three. For general use, Simpson's rule is more flexible. Like all fixed-step Newton-Cotes formulas, it is not suitable for integrands with discontinuities or singularities.
