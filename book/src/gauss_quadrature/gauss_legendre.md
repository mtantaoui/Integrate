# Gauss-Legendre

## Motivation

Instead of evaluating the integrand at equally spaced nodes as Newton-Cotes methods do,
Gauss-Legendre quadrature places nodes at the zeros of Legendre polynomials — positions
that are optimal in the sense that an \\(n\\)-point rule integrates all polynomials of
degree up to \\(2n-1\\) exactly. This means it achieves the same accuracy as a
\\(2n\\)-point trapezoidal rule with only half the function evaluations, making it the
preferred method for smooth integrands on finite intervals.

## Example

```rust,editable
use integrate::gauss_quadrature::legendre_rule;

let square = |x: f64| x * x;

let a = 0.0;
let b = 1.0;

let n: usize = 100;

let integral = legendre_rule(square, a, b, n);
println!("{}", integral);
```

## Understanding Gauss-Legendre rule

### Legendre polynomials

With respect to the inner product

\\[
\langle f, g \rangle = \int\_{-1}^{1} f(x) \cdot g(x) \cdot w(x) \\, dx
\\]

the Legendre polynomials, defined by Rodrigues' formula

\\[
P\_n(x) = \frac{1}{2^n n!} \frac{\partial^{n} (x^2 - 1)^n}{\partial x^n} \quad \text{for} \quad n > 0
\\]

and \\(P\_0(x) = 1\\), form an orthogonal family with weight function \\(w(x) = 1\\) on
\\([-1, 1]\\).

### Quadrature formula on \\([-1, 1]\\)

The \\(n\\)-point Gauss-Legendre quadrature formula \\(GL\_n(f)\\), for approximating
\\(\int\_{-1}^{1} f(x)\\,dx\\), is given by

\\[
GL\_n(f) = A\_1 f(x\_1) + \cdots + A\_n f(x\_n)
\\]

where \\(x\_1, \ldots, x\_n\\) are the zeros of \\(P\_n\\) and the weights are

\\[
A\_i = \frac{2}{1 - x\_i^2} \cdot \frac{1}{n^2 P\_{n-1}(x\_i)^2} \quad \text{for} \quad i = 1, \ldots, n
\\]

All weights \\(A\_i\\) are positive, which makes the rule numerically stable.

### Change of variables to \\([a, b]\\)

To integrate over a general interval \\([a, b]\\), the substitution
\\(t = \dfrac{2(x - a)}{b - a} - 1\\) gives

\\[
\int\_{a}^{b} f(x) \\, dx = \frac{b - a}{2} \int\_{-1}^{1} f\left( \frac{t(b-a) + (b+a)}{2} \right) dt
\\]

so the \\(n\\)-point Gauss-Legendre approximation on \\([a, b]\\) is

\\[
GL\_n(f, a, b) = A\_1' f(x\_1') + \cdots + A\_n' f(x\_n')
\\]

where the mapped nodes and weights are

\\[
x\_i' = \frac{b-a}{2} \cdot x\_i + \frac{b+a}{2}, \qquad A\_i' = \frac{b-a}{2} \cdot A\_i \quad \text{for} \quad i = 1, \ldots, n
\\]

### Truncation error

The truncation error for the rule on \\([-1, 1]\\) is

\\[
\int\_{-1}^{1} f(x) \\, dx - GL\_n(f) = K \cdot \frac{f^{(2n)}(c)}{(2n)!}
\\]

where \\(c \in (-1, 1)\\) is unknown and \\(K\\) is a constant determined by

\\[
K = \int\_{-1}^{1} x^{2n} \\, dx - GL\_n(x^{2n})
\\]

The same form holds on \\([a, b]\\) with \\(c \in (a, b)\\). A corollary is that if
\\(f^{(2n)}(x) = 0\\) for all \\(x \in [a, b]\\) — i.e. if \\(f\\) is a polynomial of
degree at most \\(2n - 1\\) — then the rule is exact.

### Node and weight computation

Most Gauss quadrature rules compute nodes and weights by building a tridiagonal Jacobi
matrix and finding its eigenvalues (Golub-Welsch algorithm), which costs \\(O(n^2)\\) in
total. The Legendre rule in this crate uses a different strategy that delivers each node
and weight independently in \\(O(1)\\) time, regardless of \\(n\\):

- **\\(n \leq 100\\):** each \\((x\_i, A\_i)\\) pair is read directly from a precomputed
  table of \\(\theta\\)-coordinates and weights stored in the source. The lookup is a
  single array index with no arithmetic at all.

- **\\(n > 100\\):** Bogaert's asymptotic formula expresses each node as
  \\(x\_i = \cos\theta\_i\\), where \\(\theta\_i\\) is obtained from the \\(k\\)-th zero of
  the Bessel function \\(J\_0\\) via a fixed-degree Chebyshev polynomial correction
  (6 terms for the node, 8 terms for the weight). The entire computation is a handful
  of fused-multiply-add operations — no Newton iteration, no eigenvalue solver.

The result is that the per-node cost is \\(O(1)\\) and the error on every node and weight
is within a few ulps of the true value, even for very large \\(n\\).

## Limitations

Gauss-Legendre quadrature is best suited for smooth functions on finite intervals. It is
not appropriate for functions with singularities, discontinuities, or functions defined on
unbounded domains. Rapidly oscillating functions may also require a very large number of
nodes \\(n\\). For semi-infinite or doubly-infinite domains use the Gauss-Laguerre or
Gauss-Hermite rules instead.

## References

Node and weight computation based on the original C++ implementation by **Ignace Bogaert**,
featuring \\(O(1)\\) per-node complexity and errors within a few ulps:

\\[
[1]: \text{Bogaert I., }
\textit{Iteration-free computation of Gauss-Legendre quadrature nodes and weights,} \\\\
\text{SIAM Journal on Scientific Computing, Volume 36, Number 3, 2014, pages A1008\text{–}A1026.}
\\]

For more details, here is a [link](https://www.cfm.brown.edu/faculty/gk/APMA2560/Handouts/GL_quad_Bogaert_2014.pdf) to the article.
