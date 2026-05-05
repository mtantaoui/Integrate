# Gauss-Hermite

## Motivation

Gauss-Hermite quadrature is designed for integrals of the form
\\(\int\_{-\infty}^{+\infty} f(x)\\,e^{-x^2}\\,dx\\), which appear naturally in quantum
mechanics, probability theory (Gaussian distributions), and signal processing. By
exploiting the orthogonality of Hermite polynomials, it provides exponential convergence
for smooth integrands with Gaussian-type decay.

## Example

```rust,editable
use integrate::gauss_quadrature::gauss_hermite_rule;

let f = |x: f64| 1.0;

let n: usize = 100;

let integral = gauss_hermite_rule(f, n);
println!("{}", integral);
```

## Understanding Gauss-Hermite rule

### Hermite polynomials

With respect to the inner product

\\[
\langle f, g \rangle = \int\_{-\infty}^{+\infty} f(x) \cdot g(x) \cdot w(x) \\, dx
\\]

the Hermite polynomials, defined by Rodrigues' formula

\\[
H\_n(x) = (-1)^n \cdot e^{x^2} \cdot \frac{\partial^{n} e^{-x^2}}{\partial x^n} \quad \text{for} \quad n > 0
\\]

and \\(H\_0(x) = 1\\), form an orthogonal family with weight function \\(w(x) = e^{-x^2}\\)
on the entire real line, satisfying the three-term recurrence

\\[
H\_{k+1}(x) = 2x\\,H\_k(x) - 2k\\,H\_{k-1}(x)
\\]

Because \\(H\_n\\) is either even or odd, its zeros are symmetric about the origin: if
\\(x\\) is a zero then so is \\(-x\\), and for odd \\(n\\) there is always a zero at
\\(x = 0\\).

### Quadrature formula

The \\(n\\)-point Gauss-Hermite formula \\(GH\_n(f)\\), approximating
\\(\int\_{-\infty}^{+\infty} f(x)\\,e^{-x^2}\\,dx\\), is

\\[
GH\_n(f) = A\_1 f(x\_1) + \cdots + A\_n f(x\_n)
\\]

where \\(x\_1, \ldots, x\_n\\) are the zeros of \\(H\_n\\) and the weights are

\\[
A\_i = \frac{2^{n+1} \cdot n! \cdot \sqrt{\pi}}{H\_{n-1}(x\_i)^2} \quad \text{for} \quad i = 1, \ldots, n
\\]

### Truncation error

The truncation error for the \\(n\\)-point rule is

\\[
\begin{split}
\int\_{-\infty}^{+\infty} f(x)\\,e^{-x^2}\\,dx - GH\_n(f) \\\\
= \frac{n!\\,\sqrt{\pi}}{2^n\\,(2n)!}\\, f^{(2n)}(\xi)
\end{split}
\\]

for some \\(\xi \in \mathbb{R}\\). A corollary is that if \\(f\\) is a polynomial of
degree at most \\(2n - 1\\) then the rule is exact.

### Node and weight computation

Nodes and weights are computed via the **Golub-Welsch algorithm**: the \\(n\\) zeros of
\\(H\_n\\) are the eigenvalues of the symmetric tridiagonal Jacobi matrix with all-zero
diagonal and off-diagonal entries

\\[
e\_k = \sqrt{\frac{k}{2}}, \quad k = 1, \ldots, n-1
\\]

The all-zero diagonal directly reflects the symmetry of \\(H\_n\\): eigenvalues come in
\\(\pm\\) pairs (with \\(0\\) as an unpaired eigenvalue for odd \\(n\\)), halving the
effective search space. Eigenvalues are located by bisection on Sturm sequences, costing
\\(O(n)\\) per node and \\(O(n^2)\\) in total.

Once each zero \\(x\_i\\) is known, \\(H\_{n-1}(x\_i)\\) is evaluated via the three-term
recurrence and the weight \\(A\_i\\) computed directly. For large \\(n\\) the factor
\\(2^{n+1} \cdot n!\\) overflows `f64`, so the implementation switches to `BigInt`
arithmetic for the numerator and denominator before converting back. A warning is printed
to `stderr` if any weight underflows to zero after the conversion.

## Limitations

Gauss-Hermite quadrature is only appropriate for integrals of the form
\\(\int\_{-\infty}^{+\infty} f(x)\\,e^{-x^2}\\,dx\\). If the integrand does not carry the
Gaussian factor \\(e^{-x^2}\\), the method will produce inaccurate results. It is not
suitable for integrands on finite intervals or those with non-Gaussian tail behavior. For
large \\(n\\) the outermost weights are extremely small and may underflow, effectively
reducing the usable order of the rule.
