# Simpson's rule

## Motivation

Simpson's rule fits a quadratic polynomial through three equally-spaced points per
subinterval and integrates it exactly, yielding fourth-order convergence for smooth
functions. It is one of the most widely used Newton-Cotes formulas and is exact for
polynomials of degree three or less, despite using only second-degree interpolants —
a consequence of the symmetric placement of nodes.

## Example

```rust,editable
use integrate::newton_cotes::simpson_rule;

let square = |x: f64| x * x;

let a = 0.0;
let b = 1.0;

let num_steps: usize = 1_000_000;

let integral = simpson_rule(square, a, b, num_steps);
println!("{}", integral);
```

## Understanding Simpson's rule

### Composite formula

Simpson's rule approximates the integral of \\(f(x)\\) on \\([a, a+h]\\) by integrating
the quadratic through \\((a, f(a))\\), \\((a+h/2, f(a+h/2))\\), and \\((a+h, f(a+h))\\).

The composite rule decomposes \\([a, b]\\) into \\(n\\) subintervals of equal length
\\(h = (b-a)/n\\) and sums over each subinterval:

\\[
\begin{align}
S\_h(f) &= \frac{h}{6} \left[ f(a) + 4f\left(a+\tfrac{h}{2}\right) + 2f(a+h) \right. \\\\
& \left. + \cdots + 2f(b-h) + 4f\left(b-\tfrac{h}{2}\right) + f(b) \right]
\end{align}
\\]

### Euler-Maclaurin expansion

An immediate consequence of the Euler-Maclaurin summation formula relates \\(S\_h(f)\\) to
the exact integral:

\\[
\begin{align}
S\_h(f) &= \int\_{a}^{b} f(x)\\,dx + \frac{h^4}{2880} \left[ f^{(3)}(b) - f^{(3)}(a) \right] - \frac{h^6}{96768} \left[ f^{(5)}(b) - f^{(5)}(a) \right] \\\\
&+ \cdots + K h^{2p-2} \left[ f^{(2p-3)}(b) - f^{(2p-3)}(a) \right] + O(h^{2p})
\end{align}
\\]

The \\(O(h^{2p})\\) remainder means that for an infinitely differentiable function whose
first \\(2p-3\\) derivatives vanish at both endpoints, the theorem guarantees only that

\\[
\lim_{h \to 0} \left| \frac{S_h(f) - \int_{a}^{b} f(x)\,dx}{h^{2p}} \right| < M
\\]

for some \\(M > 0\\), not that the rule is exact.

The expansion also shows that \\(n\\) should be chosen so that \\(h = (b-a)/n < 1\\).
For example, if \\(h = 0.1\\):

\\[
\begin{align}
S\_{0.1}(f) &= \int\_{a}^{b} f(x)\\,dx + 3.5 \cdot 10^{-8} \left[ f^{(3)}(b) - f^{(3)}(a) \right] \\\\
&- 1.033 \cdot 10^{-11} \left[ f^{(5)}(b) - f^{(5)}(a) \right] + \cdots
\end{align}
\\]

if \\(h = 0.01\\):

\\[
\begin{align}
S\_{0.01}(f) &= \int\_{a}^{b} f(x)\\,dx + 3.5 \cdot 10^{-12} \left[ f^{(3)}(b) - f^{(3)}(a) \right] \\\\
&- 1.033 \cdot 10^{-17} \left[ f^{(5)}(b) - f^{(5)}(a) \right] + \cdots
\end{align}
\\]

and if \\(h = 10\\):

\\[
\begin{align}
S\_{10}(f) &= \int\_{a}^{b} f(x)\\,dx + 3.47 \left[ f^{(3)}(b) - f^{(3)}(a) \right] \\\\
&- 10.33 \left[ f^{(5)}(b) - f^{(5)}(a) \right] + \cdots
\end{align}
\\]

### Truncation error

If \\(f \in C^4[a, b]\\), applying the mean-value theorem to the leading term of the
Euler-Maclaurin expansion gives the standard truncation error:

\\[
S\_h(f) - \int\_{a}^{b} f(x)\\,dx = \frac{h^4}{2880}(b-a)\\,f^{(4)}(c)
\\]

for some \\(c \in [a, b]\\). A corollary is that if \\(f^{(4)}(x) = 0\\) on \\([a, b]\\)
— i.e. if \\(f\\) is a polynomial of degree at most 3 — then the rule is exact. If \\(f\\)
is a cubic, \\(n\\) may be chosen to be 1.

### Computation

The \\(2n+1\\) nodes are the endpoints, the \\(n\\) subinterval midpoints, and the \\(n-1\\)
interior shared endpoints:

\\[
x\_i = a + \tfrac{i}{2}\\,h, \quad i = 0, 1, \ldots, 2n, \quad h = \frac{b-a}{n}
\\]

The effective weights are:

- \\(h/6\\) at the endpoints \\(x\_0 = a\\) and \\(x\_{2n} = b\\)
- \\(2h/3\\) at each midpoint \\(x\_1, x\_3, \ldots, x\_{2n-1}\\)
- \\(h/3\\) at each interior shared node \\(x\_2, x\_4, \ldots, x\_{2n-2}\\)

The implementation makes a single forward pass over the \\(n\\) subintervals, accumulating
three terms per step — \\(f(a) + 4f(a+h/2)\\) for the first, \\(2f(a+ih) + 4f(a+ih+h/2)\\)
for interior subintervals, and \\(f(b)\\) for the last — then multiplies the total by
\\(h/6\\). A total of \\(2n+1\\) function evaluations are required.

## Limitations

Simpson's rule converges as \\(O(h^4)\\) and requires at least four-times-differentiable
integrands for its error estimate to apply. It is not suitable for functions with
discontinuities or singularities without splitting the interval. The rule evaluates
\\(2n+1\\) points per \\(n\\) subintervals — slightly more than the trapezoidal rule's
\\(n+1\\) — so Newton's 3/8 rule may be preferred when the number of subintervals is a
multiple of three and an equivalent accuracy is sought with more uniform node spacing.
