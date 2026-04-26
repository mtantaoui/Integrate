# Trapezoidal rule

## Motivation

The trapezoidal rule approximates a definite integral by replacing the integrand with
straight-line segments between adjacent sample points. It is a practical first choice for
smooth functions and is the natural building block for Richardson extrapolation methods
such as Romberg's method.

## Example

```rust,editable
use integrate::newton_cotes::trapezoidal_rule;

let square = |x: f64| x * x;

let a = 0.0;
let b = 1.0;

let num_steps: usize = 1_000_000;

let integral = trapezoidal_rule(square, a, b, num_steps);
println!("{}", integral);
```

## Understanding the Trapezoidal rule

### Composite formula

The trapezoidal rule approximates the integral of \\(f(x)\\) on \\([a, a+h]\\) by the
signed area of the trapezoid formed by the line connecting \\((a, f(a))\\) to
\\((a+h, f(a+h))\\).

The composite rule decomposes \\([a, b]\\) into \\(n\\) subintervals of equal length
\\(h = (b-a)/n\\) and sums the result over each subinterval:

\\[
T_h(f) = h \left[ \frac{f(a)}{2} + f(a+h) + f(a+2h) + \cdots + f(b-h) + \frac{f(b)}{2} \right]
\\]

### Euler-Maclaurin expansion

The Euler-Maclaurin summation formula relates \\(T\_h(f)\\) to the exact integral:

\\[
\begin{align}
T\_h(f) &= \int\_{a}^{b} f(x)\\,dx + \frac{h^2}{12} \left[ f^\prime(b) - f^\prime(a) \right] - \frac{h^4}{720} \left[ f^{(3)}(b) - f^{(3)}(a) \right] \\\\
&+ \cdots + K h^{2p-2} \left[ f^{(2p-3)}(b) - f^{(2p-3)}(a) \right] + O(h^{2p})
\end{align}
\\]

The \\(O(h^{2p})\\) remainder means that for an infinitely differentiable function whose
first \\(2p-3\\) derivatives vanish at both endpoints, the theorem guarantees only that

\\[
\lim_{h \to 0} \left| \frac{T_h(f) - \int_{a}^{b} f(x)\,dx}{h^{2p}} \right| < M
\\]

for some \\(M > 0\\), not that the rule is exact.

The expansion also shows that \\(n\\) should be chosen so that \\(h = (b-a)/n < 1\\).
For example, if \\(h = 0.1\\):

\\[
\begin{split}
T\_{0.1}(f) &= \int\_{a}^{b} f(x)\\,dx + 0.00083 \left[ f^\prime(b) - f^\prime(a) \right] \\\\
&- 0.00000014 \left[ f^{(3)}(b) - f^{(3)}(a) \right] + \cdots
\end{split}
\\]

if \\(h = 0.01\\):

\\[
\begin{split}
T\_{0.01}(f) &= \int\_{a}^{b} f(x)\\,dx + 0.0000083 \left[ f^\prime(b) - f^\prime(a) \right] \\\\
&- 0.000000000014 \left[ f^{(3)}(b) - f^{(3)}(a) \right] + \cdots
\end{split}
\\]

and if \\(h = 10\\):

\\[
\begin{split}
T\_{10}(f) &= \int\_{a}^{b} f(x)\\,dx + 8.3333 \left[ f^\prime(b) - f^\prime(a) \right] \\\\
&- 13.89 \left[ f^{(3)}(b) - f^{(3)}(a) \right] + \cdots
\end{split}
\\]

### Truncation error

If \\(f \in C^2[a, b]\\), applying the mean-value theorem to the leading term of the
Euler-Maclaurin expansion gives the standard truncation error:

\\[
T\_h(f) - \int\_{a}^{b} f(x)\\,dx = \frac{h^2}{12}(b-a)\\,f^{\prime\prime}(c)
\\]

for some \\(c \in [a, b]\\). The positive sign reflects that the trapezoid overestimates
the integral when \\(f\\) is convex (\\(f^{\prime\prime} > 0\\)) and underestimates when
\\(f\\) is concave (\\(f^{\prime\prime} < 0\\)).

A corollary is that if \\(f^{\prime\prime}(x) = 0\\) on \\([a, b]\\) — i.e. if \\(f\\)
is linear — then the rule is exact. If \\(f\\) is linear, \\(n\\) may be chosen to be 1.

### Computation

The \\(n+1\\) nodes are the equally spaced points

\\[
x\_i = a + i\,h, \quad i = 0, 1, \ldots, n, \quad h = \frac{b-a}{n}
\\]

The endpoint nodes carry weight \\(h/2\\) and all interior nodes carry weight \\(h\\), so
the formula can be reorganised as

\\[
T\_h(f) = h\,\left(\frac{f(a) + f(b)}{2} + \sum\_{i=1}^{n-1} f(x\_i)\right)
\\]

The implementation accumulates the interior sum in a single forward pass over
\\(i = 1, \ldots, n-1\\), adds the half-weighted endpoints, then multiplies by \\(h\\).
A total of \\(n+1\\) function evaluations are required.

## Limitations

The trapezoidal rule converges as \\(O(h^2)\\) for general smooth functions. While more
accurate than the rectangle rule, it is outperformed by Simpson's rule for smooth
integrands. It should not be used for functions with discontinuities or singularities in
\\([a, b]\\) without splitting the interval first. For smooth periodic functions on
\\([a, b]\\) the rule achieves spectral convergence (faster than any power of \\(h\\))
because the endpoint correction terms in the Euler-Maclaurin expansion all vanish.
