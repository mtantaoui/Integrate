# Gauss Hermite

Gauss-Hermite quadrature formulas are used to integrate functions \\(f(x) e^{x^2}\\) from \\(-\infty\\) to \\(+\infty\\).

With respect to the inner product

\\[
\langle f,g \rangle = \int_{-\infty}^{+\infty} f(x) * g(x) * w(x) dx
\\]

the Hermite polynomials

\\[
H_n(x) = (-1)^n * e^{x^2}* \frac{\partial^{n} e^{-x^2}}{\partial x^n} \quad \text{for} \quad n > 0
\\]

and \\(H_0(x) = 1\\) form an orthogonal family of polynomials with weight function \\(w(x) = e^{-x^2}\\) on the entire \\(x\\)-axis.

The \\(n\\)-point Gauss-Hermite quadrature formula, \\(GH_n ( f(x) )\\), for approximating the integral of \\(f(x) e^{-x^2}\\) over the entire \\(x\\)-axis, is given by

\\[
GH_n ( f(x) ) = A_1 f(x_1) + ··· + A_n f(x_n)
\\]

where \\(x_i\\) , for \\(i = 1,...,n\\), are the zeros of \\(H_n\\) and

\\[
A_i = \frac{2^{n+1} * n! * \sqrt{\pi}}{H_{n-1} (x_i)^2} \quad \text{for} \quad i = 1,...,n
\\]
