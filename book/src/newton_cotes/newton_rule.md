# Newton's 3/8 rule

## Motivation

Newton's 3/8 rule fits a cubic polynomial through four equally-spaced points per
subinterval and integrates it exactly. It achieves the same order of accuracy as
Simpson's rule (\\(O(h^4)\\)) but uses three equally-spaced subdivisions per panel,
making it useful when the number of subintervals is a multiple of three or when more
uniform node spacing is preferred.

## Example

```rust,editable
use integrate::newton_cotes::newton_rule;

let square = |x: f64| x * x;

let a = 0.0;
let b = 1.0;

let num_steps: usize = 1_000_000;

let integral = newton_rule(square, a, b, num_steps);
println!("{}", integral);
```

## Understanding Newton's 3/8 rule

### Composite formula

Newton's \\(3/8\\) rule approximates the integral of \\(f(x)\\) on \\([a, a+h]\\) by
integrating the cubic through \\((a, f(a))\\), \\((a+h/3, f(a+h/3))\\),
\\((a+2h/3, f(a+2h/3))\\), and \\((a+h, f(a+h))\\).

The composite rule decomposes \\([a, b]\\) into \\(n\\) subintervals of equal length
\\(h = (b-a)/n\\) and sums over each subinterval:

\\[
\begin{align}
N\_h(f) &= \frac{h}{8} \left[ f(a) + 3f\left(a+\tfrac{h}{3}\right) + 3f\left(a+\tfrac{2h}{3}\right) + 2f(a+h) \right. \\\\
& \left. + \cdots + 2f(b-h) + 3f\left(b-\tfrac{2h}{3}\right) + 3f\left(b-\tfrac{h}{3}\right) + f(b) \right]
\end{align}
\\]

### Euler-Maclaurin expansion

An immediate consequence of the Euler-Maclaurin summation formula relates \\(N\_h(f)\\) to
the exact integral:

\\[
\begin{align}
N\_h(f) &= \int\_{a}^{b} f(x)\\,dx + \frac{h^4}{6480} \left[ f^{(3)}(b) - f^{(3)}(a) \right] - \frac{h^6}{244944} \left[ f^{(5)}(b) - f^{(5)}(a) \right] \\\\
&+ \cdots + K h^{2p-2} \left[ f^{(2p-3)}(b) - f^{(2p-3)}(a) \right] + O(h^{2p})
\end{align}
\\]

The \\(O(h^{2p})\\) remainder means that for an infinitely differentiable function whose
first \\(2p-3\\) derivatives vanish at both endpoints, the theorem guarantees only that

\\[
\lim_{h \to 0} \left| \frac{N_h(f) - \int_{a}^{b} f(x)\,dx}{h^{2p}} \right| < M
\\]

for some \\(M > 0\\), not that the rule is exact.

The expansion also shows that \\(n\\) should be chosen so that \\(h = (b-a)/n < 1\\).
For example, if \\(h = 0.1\\):

\\[
\begin{align}
N\_{0.1}(f) &= \int\_{a}^{b} f(x)\\,dx + 1.5 \cdot 10^{-8} \left[ f^{(3)}(b) - f^{(3)}(a) \right] \\\\
&- 4.1 \cdot 10^{-12} \left[ f^{(5)}(b) - f^{(5)}(a) \right] + \cdots
\end{align}
\\]

if \\(h = 0.01\\):

\\[
\begin{align}
N\_{0.01}(f) &= \int\_{a}^{b} f(x)\\,dx + 1.5 \cdot 10^{-12} \left[ f^{(3)}(b) - f^{(3)}(a) \right] \\\\
&- 4.1 \cdot 10^{-18} \left[ f^{(5)}(b) - f^{(5)}(a) \right] + \cdots
\end{align}
\\]

and if \\(h = 10\\):

\\[
\begin{align}
N\_{10}(f) &= \int\_{a}^{b} f(x)\\,dx + 1.54 \left[ f^{(3)}(b) - f^{(3)}(a) \right] \\\\
&- 4.08 \left[ f^{(5)}(b) - f^{(5)}(a) \right] + \cdots
\end{align}
\\]

### Truncation error

If \\(f \in C^4[a, b]\\), applying the mean-value theorem to the leading term of the
Euler-Maclaurin expansion gives the standard truncation error:

\\[
N\_h(f) - \int\_{a}^{b} f(x)\\,dx = \frac{h^4}{6480}(b-a)\\,f^{(4)}(c)
\\]

for some \\(c \in [a, b]\\). A corollary is that if \\(f^{(4)}(x) = 0\\) on \\([a, b]\\)
— i.e. if \\(f\\) is a polynomial of degree at most 3 — then the rule is exact. If \\(f\\)
is a cubic, \\(n\\) may be chosen to be 1.

### Computation

The \\(3n+1\\) nodes are the endpoints and the two third-points of each subinterval:

\\[
x\_{3i+k} = a + \left(i + \tfrac{k}{3}\right)h, \quad
i = 0, \ldots, n-1,\quad k = 0, 1, 2, \quad
\text{plus } x\_{3n} = b
\\]

The effective weights are:

- \\(h/8\\) at the endpoints \\(x\_0 = a\\) and \\(x\_{3n} = b\\)
- \\(3h/8\\) at each third-point \\(x\_1, x\_2, x\_4, x\_5, \ldots\\)
- \\(h/4\\) at each interior shared node \\(x\_3, x\_6, \ldots, x\_{3n-3}\\)

The implementation makes a single forward pass over the \\(n\\) subintervals, accumulating
four terms per step — \\(f(a) + 3f(a+h/3) + 3f(a+2h/3)\\) for the first, \\(2f(a+ih) +
3f(a+ih+h/3) + 3f(a+ih+2h/3)\\) for interior subintervals, and \\(f(b)\\) for the last —
then multiplies the total by \\(h/8\\). A total of \\(3n+1\\) function evaluations are
required.

## Limitations

Newton's 3/8 rule converges as \\(O(h^4)\\), the same order as Simpson's rule, but
requires \\(3n+1\\) evaluations versus Simpson's \\(2n+1\\) for the same \\(n\\). For
general use, Simpson's rule is more efficient. Like all fixed-step Newton-Cotes formulas,
it is not suitable for integrands with discontinuities or singularities inside \\([a, b]\\)
without splitting the interval first.
