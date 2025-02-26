# Gauss-Legendre

Gauss-Legendre quadrature formulas are used to integrate functions \\(f(x)\\) over a closed and bounded interval \\(\[a, b\]\\).

Let \\(\int\_{a}^{b} f(x) dx\\) denote the integral of \\(f(x)\\) from \\(a\\) to \\(b\\). After making the change of variable \\(t = \dfrac{2(x-a)}{b-a} - 1\\), then

\\[
\int_{a}^{b} f(x) dx = \frac{b-a}{2} \int_{-1}^{1} f \left( \frac{t(b-a) + (b+a)}{2} \right) dt
\\]

with respect to the inner product

\\[
\langle f,g \rangle = \int\_{-1}^{1} f(x) \cdot g(x) \cdot w(x) dx
\\]

The Legendre polynomials are defined by

\\[
P_n(x) = \frac{1}{2^n n! } \frac{\partial^{n} (x^2 - 1)^n}{\partial x^n} \quad \text{for} \quad n>0
\\]

and \\(P_0(x) = 1\\), form an orthogonal family of polynomials with weight function \\(w(x) = 1\\) on the interval \\([-1,1]\\).

The \\(n\\)-point Gauss-Legendre quadrature formula, \\(GL_n ( f )\\), for approximating the integral of \\(f(x)\\) over \\([-1,1]\\), is given by

\\[
GL_n ( f ) = A_1 f(x_1) + \cdots + A_n f(x_n)
\\]

where \\(x_i\\), \\(i = 1,\dots,n\\), are the zeros of \\(P_n\\) and

\\[
A_i = 2 \cdot \frac{1 - x_i^2}{n^2 P\_{n-1} (x_i)^2} \quad \text{for} \quad i = 1,\dots,n
\\]

The truncation error is

\\[
\int\_{-1}^{1} f(x) dx - GL_n(f) = K \cdot \frac{ f^{(2n)}(c) }{2n!}
\\]

where \\(K\\) is a constant, and \\(c\\) is some unknown number that verifies \\(-1 < c < 1\\). The constant \\(K\\) is easily determined from

\\[
K = \int\_{-1}^{1} x^{2n} dx - GL_n (x^{2n} )
\\]

Generalizing, in order to integrate \\(f(x)\\) over \\([a,b]\\), the \\(n\\)-point Gauss-Legendre quadrature formula, \\(GL_n ( f(x), a, b )\\), is given by

\\[
GL_n ( f(x), a, b ) = A_1' f(x_1') + \cdots + A_n' f(x_n') \quad \text{where} \quad x_i' = \frac{b-a}{2} \cdot x_i + \frac{b+a}{2}
\\]

\\(x_i\\), \\(i = 1,\dots,n\\), are the zeros of \\(P_n\\) and

\\[
A_i' = (b-a) \frac{1 - x_i^2}{n^2 P\_{n-1}(x_i)^2} = \frac{b-a}{2} \cdot A_i, \quad i = 1,\dots,n
\\]

The truncation error is

\\[
\int\_{a}^{b} f(x) dx - GL_n(f, a, b) = K \cdot \frac{ f^{(2n)}(c) }{2n!}
\\]

where \\(K\\) is a constant, and \\(c\\) is some unknown number that verifies \\(a < c < b\\).

## References:

Original C++ code implemented by **Ignace Bogaert**.

The main features of this software are:

- **Speed**: due to the simple formulas and the \\(O(1)\\) complexity computation of
  individual Gauss-Legendre quadrature nodes and weights.
  This makes it compatible with parallel computing paradigms.
- **Accuracy**: the error on the nodes and weights is within a few ulps.

### Ignace Bogaert's Paper:

\\[
[1]: \text{Ignace Bogaert,}
\textit{ Iteration-free computation of Gauss-Legendre quadrature nodes and weights,} \\\\
\text{ SIAM Journal on Scientific Computing, Volume 36, Number 3, 2014, pages A1008-1026.}
\\]

For more details, here is a [link](https://www.cfm.brown.edu/faculty/gk/APMA2560/Handouts/GL_quad_Bogaert_2014.pdf) to the article.
