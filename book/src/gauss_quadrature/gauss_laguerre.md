# Gauss-Laguerre

## Motivation

Gauss-Laguerre quadrature is designed for integrals of the form
\\(\int\_0^{+\infty} f(x)\\,e^{-x}\\,dx\\), where the exponential decay naturally arises in
probability distributions, physics, and engineering problems. By using the zeros and
weights of Laguerre polynomials, the method achieves high accuracy with far fewer
function evaluations than general-purpose rules on a truncated domain.

## Example

```rust,editable
use integrate::gauss_quadrature::gauss_laguerre_rule;

let f = |x: f64| 1.0;

let n: usize = 100;

let integral = gauss_laguerre_rule(f, n);
println!("{}", integral);
```

## Understanding Gauss-Laguerre rule

### Laguerre polynomials

With respect to the inner product

\\[
\langle f, g \rangle = \int\_{0}^{+\infty} f(x) \cdot g(x) \cdot w(x) \\, dx
\\]

the Laguerre polynomials, defined by Rodrigues' formula

\\[
L\_n(x) = e^x \frac{\partial^{n} (x^n e^{-x})}{\partial x^n} \quad \text{for} \quad n > 0
\\]

and \\(L\_0(x) = 1\\), form an orthogonal family with weight function \\(w(x) = e^{-x}\\)
on the positive \\(x\\)-axis, satisfying the three-term recurrence

\\[
L\_{k+1}(x) = \frac{(2k + 1 - x)\\,L\_k(x) - k\\,L\_{k-1}(x)}{k + 1}
\\]

### Quadrature formula

The \\(n\\)-point Gauss-Laguerre formula \\(GL\_n(f)\\), approximating
\\(\int\_0^{+\infty} f(x)\\,e^{-x}\\,dx\\), is

\\[
GL\_n(f) = A\_1 f(x\_1) + \cdots + A\_n f(x\_n)
\\]

where \\(x\_1, \ldots, x\_n\\) are the zeros of \\(L\_n\\) and the weights are

\\[
A\_i = \frac{x\_i}{(n+1)^2 \\, L\_{n+1}(x\_i)^2} \quad \text{for} \quad i = 1, \ldots, n
\\]

### Truncation error

The truncation error for the \\(n\\)-point rule is

\\[
\int\_0^{+\infty} f(x)\\,e^{-x}\\,dx - GL\_n(f) = \frac{(n!)^2}{(2n)!}\\, f^{(2n)}(\xi)
\\]

for some \\(\xi > 0\\). A corollary is that if \\(f\\) is a polynomial of degree at most
\\(2n - 1\\) then the rule is exact.

### Node and weight computation

Nodes and weights are computed via the **Golub-Welsch algorithm**: the \\(n\\) zeros of
\\(L\_n\\) are the eigenvalues of the symmetric tridiagonal Jacobi matrix

\\[
J = \begin{pmatrix}
d\_0 & e\_1 & & \\\\
e\_1 & d\_1 & e\_2 & \\\\
& \ddots & \ddots & e\_{n-1} \\\\
& & e\_{n-1} & d\_{n-1}
\end{pmatrix}
\\]

with diagonal entries \\(d\_k = 2k + 1\\) and off-diagonal entries \\(e\_k = k\\) for
\\(k = 0, 1, \ldots, n-1\\).

Eigenvalues are located by bisection on Sturm sequences, costing \\(O(n)\\) per node and
\\(O(n^2)\\) in total. Once each zero \\(x\_i\\) is known, the three-term recurrence is
used to evaluate \\(L\_{n+1}(x\_i)\\) and compute the weight \\(A\_i\\) directly.

For large \\(n\\) the weights of the outermost nodes decrease extremely rapidly and may
underflow to zero; a warning is printed to `stderr` when this occurs.

## Limitations

Gauss-Laguerre quadrature is only appropriate for integrals of the form
\\(\int\_0^{+\infty} f(x)\\,e^{-x}\\,dx\\). If the integrand does not carry the
exponential factor \\(e^{-x}\\), the method will produce inaccurate results. It is not
suitable for integrands defined on finite intervals or on the full real line.
