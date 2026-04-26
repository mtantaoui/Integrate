# Rectangle rule

## Motivation

The rectangle rule (also called the midpoint rule) is the simplest Newton-Cotes formula
for numerical integration. It approximates each subinterval's contribution by the
function value at its midpoint, making it easy to implement and exact for constant
functions. Despite its low order of accuracy, it is a useful baseline and requires no
endpoint evaluations.

## Example

```rust,editable
use integrate::newton_cotes::rectangle_rule;

let square = |x: f64| x * x;

let a = 0.0;
let b = 1.0;

let num_steps: usize = 1_000_000;

let integral = rectangle_rule(square, a, b, num_steps);
println!("{}", integral);
```

## Understanding the Rectangle rule

### Composite formula

The rectangle rule approximates the integral of \\(f(x)\\) on \\([a, a+h]\\) by the
signed area of the rectangle with width \\(h\\) and height \\(f\left(a + h/2\right)\\).

The composite rule decomposes \\([a, b]\\) into \\(n\\) subintervals of equal length
\\(h = (b-a)/n\\) and sums the result over each subinterval:

\\[
R\_h(f) = h \left[ f\left(a + \tfrac{h}{2}\right) + f\left(a + \tfrac{3h}{2}\right) + \cdots + f\left(b - \tfrac{h}{2}\right) \right]
\\]

### Euler-Maclaurin expansion

The Euler-Maclaurin summation formula relates \\(R\_h(f)\\) to the exact integral:

\\[
\begin{align}
R\_h(f) &= \int\_{a}^{b} f(x)\\,dx - \frac{h^2}{24} \left[ f^\prime(b) - f^\prime(a) \right] + \frac{7h^4}{5760} \left[ f^{(3)}(b) - f^{(3)}(a) \right] \\\\
&+ \cdots + K h^{2p-2} \left[ f^{(2p-3)}(b) - f^{(2p-3)}(a) \right] + O(h^{2p})
\end{align}
\\]

The \\(O(h^{2p})\\) remainder means that for an infinitely differentiable function whose
first \\(2p-3\\) derivatives vanish at both endpoints, the theorem guarantees only that

\\[
\lim_{h \to 0} \left| \frac{R_h(f) - \int_{a}^{b} f(x)\,dx}{h^{2p}} \right| < M
\\]

for some \\(M > 0\\), not that the rule is exact.

The expansion also shows that \\(n\\) should be chosen so that \\(h = (b-a)/n < 1\\).
For example, if \\(h = 0.1\\):

\\[
\begin{split}
R\_{0.1}(f) &= \int\_{a}^{b} f(x)\\,dx - 0.00042 \left[ f^\prime(b) - f^\prime(a) \right] \\\\
&+ 0.00000012 \left[ f^{(3)}(b) - f^{(3)}(a) \right] + \cdots
\end{split}
\\]

if \\(h = 0.01\\):

\\[
\begin{split}
R\_{0.01}(f) &= \int\_{a}^{b} f(x)\\,dx - 0.0000042 \left[ f^\prime(b) - f^\prime(a) \right] \\\\
&+ 0.000000000012 \left[ f^{(3)}(b) - f^{(3)}(a) \right] + \cdots
\end{split}
\\]

and if \\(h = 10\\):

\\[
\begin{split}
R\_{10}(f) &= \int\_{a}^{b} f(x)\\,dx - 4.1667 \left[ f^\prime(b) - f^\prime(a) \right] \\\\
&+ 12.15 \left[ f^{(3)}(b) - f^{(3)}(a) \right] + \cdots
\end{split}
\\]

### Truncation error

If \\(f \in C^2[a, b]\\), applying the mean-value theorem to the leading term of the
Euler-Maclaurin expansion gives the standard truncation error:

\\[
R\_h(f) - \int\_{a}^{b} f(x)\\,dx = -\frac{h^2}{24}(b-a)\\,f^{\prime\prime}(c)
\\]

for some \\(c \in [a, b]\\). A corollary is that if \\(f^{\prime\prime}(x) = 0\\) on
\\([a, b]\\) — i.e. if \\(f\\) is linear — then the rule is exact. If \\(f\\) is linear,
\\(n\\) may be chosen to be 1.

### Computation

The \\(n\\) nodes are the midpoints of the \\(n\\) subintervals:

\\[
x\_i = a + \left(i + \tfrac{1}{2}\right)h, \quad i = 0, 1, \ldots, n-1, \quad h = \frac{b-a}{n}
\\]

Each node carries weight \\(h\\), so the sum is:

\\[
R\_h(f) = h \sum\_{i=0}^{n-1} f(x\_i)
\\]

The implementation makes a single forward pass over the \\(n\\) midpoints, accumulating
the weighted sum, then multiplies by \\(h\\). No endpoint evaluations are needed.

## Limitations

The rectangle rule converges as \\(O(h^2)\\), making it the least accurate of the
Newton-Cotes rules in this library. Prefer Simpson's rule or Romberg's method for
smooth integrands. The rule is not suitable for functions with discontinuities or
singularities inside \\([a, b]\\) without splitting the interval first.
