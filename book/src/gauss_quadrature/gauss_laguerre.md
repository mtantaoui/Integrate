# Gauss-Laguerre

## Motivation

Gauss-Laguerre quadrature is designed for integrals of the form \\(\int_0^{+\infty} f(x)\\,e^{-x}\\,dx\\), where the exponential decay naturally arises in probability distributions, physics, and engineering problems. By using the zeros and weights of Laguerre polynomials, the method achieves high accuracy with far fewer function evaluations than general-purpose rules on a truncated domain.

## Example

```rust,editable
use integrate::gauss_quadrature::gauss_laguerre_rule;

let f = |x: f64| 1.0;

let n:usize = 100;

let integral = gauss_laguerre_rule(f, n);
println!("{}",integral);
```

## Understanding Gauss-Laguerre rule

Gauss-Laguerre quadrature formulas are used to integrate functions \\(f(x) e^{-x}\\) over the positive \\(x\\)-axis.

With respect to the inner product

\\[
\langle f,g \rangle = \int\_{0}^{+\infty} f(x) \cdot g(x) \cdot w(x) dx
\\]

the Laguerre polynomials are defined by

\\[
L_n(x) = e^x \dfrac{\partial^{n} x^n e^{-x}}{\partial x^n}, \quad \text{for} \quad n > 0
\\]

and \\(L_0(x) = 1\\) form an orthogonal family of polynomials with weight function \\(w(x) = e^{-x}\\) on the positive \\(x\\)-axis.

The \\(n\\)-point Gauss-Laguerre quadrature formula, \\(GL_n ( f(x) )\\), for approximating the integral of \\(f(x) e^{-x}\\) over \\(\left[0, \infty \right[\\), is given by

\\[
GL_n ( f(x) ) = A_1 f(x_1) + \cdots + A_n f(x_n)
\\]

where \\(x_i\\), \\(i = 1,\dots,n\\), are the zeros of \\(L_n\\) and

\\[
A_i = \dfrac{n!^2}{ x_i L\_{n-1} (x_i)^2} \quad \text{for} \quad i = 1,\dots,n
\\]

## Limitations

Gauss-Laguerre quadrature is only appropriate for integrals of the form \\(\int_0^{+\infty} f(x)\\,e^{-x}\\,dx\\). If the integrand does not have exponential decay, the method will produce inaccurate results. It is also not suitable for integrands defined on finite intervals or on the full real line.
