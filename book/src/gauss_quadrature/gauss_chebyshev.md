# Gauss-Chebyshev Quadrature

## Motivation

Chebyshev quadrature is designed for integrands that carry an algebraic singularity at
the endpoints of \\([-1, 1]\\). Rather than struggling with the singularity numerically,
the weight function absorbs it exactly into the orthogonality relation, leaving a smooth
\\(f(x)\\) to be sampled.

Two weight functions are supported:

- **First kind** — \\(w(x) = 1/\sqrt{1-x^2}\\), which diverges at \\(\pm 1\\) and
  arises in arclength integrals, potential theory, and spectral methods.
- **Second kind** — \\(w(x) = \sqrt{1-x^2}\\), which vanishes at \\(\pm 1\\) and
  appears in problems with elliptic geometry or Jacobi-weight boundary conditions.

Because the zeros and weights of both families have exact trigonometric closed forms,
node and weight computation requires only a single trigonometric evaluation per point
— making this the most computationally efficient Gauss rule in the library.

---

## Chebyshev Polynomials of the First Kind

### Example

```rust,editable
use integrate::gauss_quadrature::gauss_first_kind_chebyshev_rule;

let f = |x: f64| 1.0;

let n: usize = 100;

let integral = gauss_first_kind_chebyshev_rule(f, n);
println!("{}", integral);
```

### Understanding Gauss-Chebyshev rule (first kind)

With respect to the inner product

\\[
\langle f, g \rangle = \int\_{-1}^{1} f(x) \cdot g(x) \cdot w(x) \\, dx
\\]

the Chebyshev polynomials of the first kind, defined by

\\[
T\_n(x) = \cos(n \arccos(x)) \quad \text{for} \quad n > 0
\\]

and \\(T\_0(x) = 1\\), form an orthogonal family with weight function
\\(w(x) = 1/\sqrt{1-x^2}\\) on \\([-1, 1]\\), satisfying the three-term recurrence

\\[
T\_{k+1}(x) = 2x\\,T\_k(x) - T\_{k-1}(x)
\\]

The \\(n\\)-point Gauss-Chebyshev quadrature formula of the first kind, \\(GC\_n^{(1)}(f)\\),
approximating \\(\int\_{-1}^{1} \frac{f(x)}{\sqrt{1-x^2}}\\,dx\\), is

\\[
GC\_n^{(1)}(f) = A\_1 f(x\_1) + \cdots + A\_n f(x\_n)
\\]

where \\(x\_1, \ldots, x\_n\\) are the zeros of \\(T\_n\\) and all weights are equal:

\\[
A\_i = \frac{\pi}{n} \quad \text{for} \quad i = 1, \ldots, n
\\]

### Truncation error

\\[
\begin{split}
\int\_{-1}^{1} \frac{f(x)}{\sqrt{1-x^2}}\\,dx - GC\_n^{(1)}(f) \\\\
= \frac{\pi}{2^{2n-1}\\,(2n)!}\\, f^{(2n)}(\xi)
\end{split}
\\]

for some \\(\xi \in (-1, 1)\\). The rule is exact for all polynomials of degree at most
\\(2n - 1\\).

### Node and weight computation

The zeros of \\(T\_n\\) have the explicit closed form

\\[
x\_k = \cos\left(\frac{(2k-1)\,\pi}{2n}\right), \quad k = 1, \ldots, n
\\]

and all weights are the same constant \\(\pi/n\\). Each node requires a single cosine
evaluation; no weight computation is needed beyond a single division. There is no
Jacobi matrix, no eigenvalue solver, and no iterative refinement — the total cost is
\\(O(n)\\) with a very small constant.

---

## Chebyshev Polynomials of the Second Kind

### Example

```rust,editable
use integrate::gauss_quadrature::gauss_second_kind_chebyshev_rule;

let f = |x: f64| 1.0;

let n: usize = 100;

let integral = gauss_second_kind_chebyshev_rule(f, n);
println!("{}", integral);
```

### Understanding Gauss-Chebyshev rule (second kind)

With respect to the same inner product structure, the Chebyshev polynomials of the
second kind, defined by

\\[
U\_n(x) = \frac{\sin((n+1)\arccos(x))}{\sin(\arccos(x))} \quad \text{for} \quad n > 0
\\]

and \\(U\_0(x) = 1\\), form an orthogonal family with weight function
\\(w(x) = \sqrt{1-x^2}\\) on \\([-1, 1]\\), satisfying the same three-term recurrence

\\[
U\_{k+1}(x) = 2x\\,U\_k(x) - U\_{k-1}(x)
\\]

The \\(n\\)-point Gauss-Chebyshev quadrature formula of the second kind, \\(GC\_n^{(2)}(f)\\),
approximating \\(\int\_{-1}^{1} f(x)\sqrt{1-x^2}\\,dx\\), is

\\[
GC\_n^{(2)}(f) = A\_1 f(x\_1) + \cdots + A\_n f(x\_n)
\\]

where \\(x\_1, \ldots, x\_n\\) are the zeros of \\(U\_n\\) and the weights are

\\[
A\_i = \frac{\pi}{n+1} \sin^2\left(\frac{i\,\pi}{n+1}\right) \quad \text{for} \quad i = 1, \ldots, n
\\]

### Truncation error

\\[
\begin{split}
\int\_{-1}^{1} f(x)\sqrt{1-x^2}\\,dx - GC\_n^{(2)}(f) \\\\
= \frac{\pi}{2^{2n+1}\\,(2n)!}\\, f^{(2n)}(\xi)
\end{split}
\\]

for some \\(\xi \in (-1, 1)\\). The rule is exact for all polynomials of degree at most
\\(2n - 1\\).

### Node and weight computation

The zeros of \\(U\_n\\) have the explicit closed form

\\[
x\_k = \cos\left(\frac{k\,\pi}{n+1}\right), \quad k = 1, \ldots, n
\\]

and each weight requires one additional sine evaluation:

\\[
A\_k = \frac{\pi}{n+1} \sin^2\left(\frac{k\,\pi}{n+1}\right)
\\]

As with the first kind, there is no Jacobi matrix, no eigenvalue solver, and no
iterative refinement. The total cost is \\(O(n)\\) with a very small constant.

The equal-angle spacing of \\(k\pi/(n+1)\\) is the same grid used in the discrete sine
transform (DST), so the second-kind nodes and weights can also be interpreted as a
DST-I evaluation.

---

## Limitations

Gauss-Chebyshev quadrature is only suitable when the integrand can be expressed in the
weighted form \\(f(x)/\sqrt{1-x^2}\\) (first kind) or \\(f(x)\sqrt{1-x^2}\\) (second
kind) on \\([-1, 1]\\). Applying either rule to an unweighted integrand will produce
inaccurate results, since the accuracy guarantee relies on orthogonality with respect to
the corresponding weight function. Both rules are restricted to the interval \\([-1, 1]\\);
for other domains or weight functions use a different Gauss rule.
