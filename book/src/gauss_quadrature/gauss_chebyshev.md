# Gauss-Chebyshev Quadrature

Gauss-Chebyshev quadrature formulas are used to integrate functions like \\(\frac{f(x)}{\sqrt{1- x^2}}\\) and \\(f(x) \cdot \sqrt{1- x^2}\\) from \\(-1\\) to \\(1\\).

## Chebyshev Polynomials of the First Kind

With respect to the inner product

\\[
\langle f,g \rangle = \int\_{-\infty}^{+\infty} f(x) \cdot g(x) \cdot w(x) dx
\\]

the Chebyshev polynomials are defined by

\\[
T_n(x) = \cos(n \cdot \arccos(x)) \quad \text{for} \quad n>0
\\]

and \\(T_0(x)=1\\) form an orthogonal family of polynomials with weight function \\(w(x)=\frac{1}{\sqrt{1 - x^2}}\\) on \\([-1, 1]\\).

The \\(n\\)-point Gauss-Chebyshev quadrature formula, \\(GC_n(f(x))\\), for approximating the integral of \\(\frac{f(x)}{\sqrt{1 - x^2}}\\) over \\([-1, 1]\\), is given by

\\[
GC_n ( f(x) ) = A_1 f(x_1) + \cdots + A_n f(x_n)
\\]

where \\(x_i\\), \\(i = 1,\dots,n\\), are the zeros of \\(T_n\\) and \\(A_i = \frac{\pi}{n}\\), \\(i = 1,\dots,n\\).

## Chebyshev Polynomials of the Second Kind

With respect to the inner product

\\[
\langle f,g \rangle = \int\_{-\infty}^{+\infty} f(x) \cdot g(x) \cdot w(x) dx
\\]

the Chebyshev polynomials are defined by

\\[
U_n(x) \cdot \sin(\arccos(x)) = \sin((n+1) \cdot \arccos(x)) \quad \text{for} \quad n>0
\\]

and \\(U_0(x)=1\\) form an orthogonal family of polynomials with weight function \\(w(x)=\sqrt{1 - x^2}\\) on \\([-1, 1]\\).

The \\(n\\)-point Gauss-Chebyshev quadrature formula, \\(GC_n(f(x))\\), for approximating the integral of \\(f(x) \cdot \sqrt{1 - x^2}\\) over \\([-1, 1]\\), is given by

\\[
GC_n ( f(x) ) = A_1 f(x_1) + \cdots + A_n f(x_n)
\\]

where \\(x_i\\), \\(i = 1,\dots,n\\), are the zeros of \\(U_n\\) and

\\[
A_i = \frac{\pi}{n + 1} \cdot \sin^2\left(\frac{i\pi}{n + 1} \right) \quad \text{for} \quad i = 1,\dots,n.
\\]
