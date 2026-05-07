# Gauss-Kronrod

## Motivation

The \\((2N+1)\\)-point Gauss-Kronrod rule simultaneously produces a high-accuracy integral estimate
and a reliable error bound in a single pass over \\(2N+1\\) function evaluations. It does so by
augmenting an \\(N\\)-point Gauss-Legendre rule with \\(N+1\\) optimally placed Kronrod points: the
original \\(N\\) Gauss nodes are reused, so no function evaluations are wasted.

This makes it the method of choice when both accuracy and a trustworthy error estimate are needed.
It is the engine behind [QUADPACK](https://en.wikipedia.org/wiki/QUADPACK) and SciPy's `scipy.integrate.quad`.

## Example

```rust,editable
use integrate::gauss_kronrod::gauss_kronrod_rule;

let f = |x: f64| x.exp();

let a = 0.0_f64;
let b = 1.0_f64;
let n = 7_usize; // 15-point Gauss-Kronrod rule

match gauss_kronrod_rule(f, a, b, n) {
    Ok((integral, error)) => {
        println!("Integral:       {:.10}", integral);
        println!("Error estimate: {:.2e}", error);
    }
    Err(e) => println!("Error: {}", e),
}
```

## Understanding Gauss-Kronrod rule

The Gauss-Kronrod rule is constructed in two stages.

### Stage 1: Gauss-Legendre base rule

The \\(N\\)-point Gauss-Legendre rule approximates \\(\int\_{-1}^{1} f(x)\\,dx\\) by

\\[
G\_N(f) = \sum\_{i=1}^{N} w\_i^G f(x\_i)
\\]

where \\(x\_1, \ldots, x\_N\\) are the zeros of the Legendre polynomial \\(P\_N\\) and
\\(w\_1^G, \ldots, w\_N^G\\) are positive weights chosen so that the formula is exact for all
polynomials of degree at most \\(2N - 1\\).

### Stage 2: Kronrod extension

The Stieltjes polynomial \\(E\_{N+1}\\) of degree \\(N+1\\) is the unique monic polynomial satisfying

\\[
\int\_{-1}^{1} E\_{N+1}(x) \cdot x^k \cdot P\_N(x) \\, dx = 0, \quad k = 0, 1, \ldots, N
\\]

Its \\(N+1\\) zeros \\(t\_1, \ldots, t\_{N+1}\\) interlace strictly with those of \\(P\_N\\) and all lie
in \\((-1, 1)\\). Adding them to the \\(N\\) existing Gauss-Legendre nodes gives \\(2N+1\\)
distinct abscissas. Weights for all \\(2N+1\\) points can then be chosen so that the combined
formula

\\[
GK\_{2N+1}(f) = \sum\_{i=1}^{N} w\_i^{GK} f(x\_i) + \sum\_{j=1}^{N+1} v\_j^{GK} f(t\_j)
\\]

is exact for all polynomials of degree at most \\(3N + 1\\) (when \\(N\\) is odd). Since the
Gauss nodes are reused, the total cost is \\(2N+1\\) function evaluations.

### Change of variables to \\([a, b]\\)

For an integral over a general interval \\([a, b]\\), the substitution
\\(x = \dfrac{(b-a)\\,t + (b+a)}{2}\\) gives

\\[
\int\_a^b f(x) \\, dx = \frac{b-a}{2} \int\_{-1}^{1} f\left(\frac{(b-a)\\,t + (b+a)}{2}\right) dt
\\]

so the \\((2N+1)\\)-point Gauss-Kronrod approximation on \\([a, b]\\) is

\\[
GK\_{2N+1}(f, a, b) = \frac{b-a}{2} \cdot GK\_{2N+1}(\tilde{f}), \quad \tilde{f}(t) = f\left(\frac{(b-a)\\,t + (b+a)}{2}\right)
\\]

### Error estimate

`gauss_kronrod_rule` returns both the integral and an error estimate:

- **integral** \\(= GK\_{2N+1}(f, a, b)\\)
- **error** \\(= \left\lvert GK\_{2N+1}(f, a, b) - G\_N(f, a, b) \right\rvert\\)

The error is the absolute difference between the \\((2N+1)\\)-point Kronrod estimate and the
\\(N\\)-point Gauss-Legendre estimate computed on the same set of function values. For smooth
integrands this quantity is a reliable indicator of the absolute integration error.

### Choosing \\(N\\)

| \\(N\\) | Total points | Degree of exactness (\\(N\\) odd: \\(3N+1\\)) |
|---------|-------------|-----------------------------------------------|
| 7       | 15          | 22                                            |
| 10      | 21          | 31                                            |
| 15      | 31          | 46                                            |

The default \\(N = 7\\) (15-point rule) matches the default used by QUADPACK and SciPy on each
subinterval and gives double-precision accuracy for most smooth integrands.

## Limitations

The single-interval Gauss-Kronrod rule requires the integrand to be smooth on \\([a, b]\\). For
integrands with singularities, discontinuities, or sharp peaks the error estimate may be
unreliable and the accuracy may be poor. In such cases the interval should be split manually
or an adaptive method should be used instead.

The rule is defined only on finite intervals \\([a, b]\\). For semi-infinite or doubly-infinite
domains use the Gauss-Laguerre or Gauss-Hermite rules instead.

## References

Node and weight computation ported from the original Fortran 77 by **Robert Piessens** and
**Maria Branders** (C translation by **John Burkardt**, GNU LGPL):

\\[
[1]: \text{Piessens R., Branders M., }
\textit{A Note on the Optimal Addition of Abscissas to Quadrature Formulas of Gauss and Lobatto,} \\\\
\text{Mathematics of Computation, Vol. 28, No. 125, 1974, pp. 135\text{–}139.}
\\]
