# Romberg's method

## Motivation

Romberg's method combines the composite trapezoidal rule with Richardson extrapolation
to systematically cancel the leading error terms from the Euler-Maclaurin expansion.
Each extrapolation step eliminates one error term, so with \\(n\\) columns the final
estimate has accuracy \\(O(h^{2n})\\) — matching a Gauss rule of order \\(2n\\) but
using only \\(2^{n-1}+1\\) function evaluations. It is an excellent choice when high
precision is required for smooth integrands on finite intervals.

## Example

```rust,editable
use integrate::romberg::romberg_method;

fn square(x: f64) -> f64 {
    x.powi(2)
}

let a = 0.0;
let b = 1.0;

let num_steps: usize = 10;

let integral = romberg_method(square, a, b, num_steps);
println!("{}", integral);
```

## Understanding Romberg's method

### Euler-Maclaurin expansion

The Euler-Maclaurin summation formula expresses the composite trapezoidal estimate
\\(T\_h(f)\\) with step size \\(h = (b-a)/n\\) as

\\[
\begin{align}
T\_h(f) &= \int\_{a}^{b} f(x)\\,dx + \frac{h^2}{12} \left[ f^\prime(b) - f^\prime(a) \right] - \frac{h^4}{720} \left[ f^{(3)}(b) - f^{(3)}(a) \right] \\\\
&+ \cdots + K h^{2p-2} \left[ f^{(2p-3)}(b) - f^{(2p-3)}(a) \right] + O(h^{2p})
\end{align}
\\]

Two features make this expansion ideal for Richardson extrapolation:

1. The error is a **pure even-power series** in \\(h\\): terms \\(h^2, h^4, h^6, \ldots\\)
   with coefficients that depend on the endpoints but not on \\(h\\).
2. When the step size is halved to \\(h/2\\), each coefficient is divided by the
   corresponding power of 4 — which allows exact cancellation.

### Richardson extrapolation

Halving \\(h\\) gives a new estimate

\\[
\begin{split}
T\_{h/2}(f) &= \int\_{a}^{b} f(x)\\,dx + \frac{h^2}{4 \cdot 12} \left[ f^\prime(b) - f^\prime(a) \right] \\\\
&- \frac{h^4}{16 \cdot 720} \left[ f^{(3)}(b) - f^{(3)}(a) \right] + \cdots
\end{split}
\\]

Forming the linear combination \\(\frac{4 T\_{h/2}(f) - T\_h(f)}{3}\\) cancels the
\\(h^2\\) term exactly:

\\[
\begin{split}
\frac{4 T\_{h/2}(f) - T\_h(f)}{3} &= \int\_{a}^{b} f(x)\\,dx + \frac{h^4}{2880} \left[ f^{(3)}(b) - f^{(3)}(a) \right] \\\\
&- \frac{h^6}{96768} \left[ f^{(5)}(b) - f^{(5)}(a) \right] + \cdots
\end{split}
\\]

This extrapolation step improves the order from \\(O(h^2)\\) to \\(O(h^4)\\). Applying it
again to two such \\(O(h^4)\\) estimates cancels the \\(h^4\\) term and yields \\(O(h^6)\\),
and so on.

### The Romberg table

The full process is organised into a triangular table \\(R[i, j]\\), where \\(i\\) is the
**refinement level** (step size \\(h\_i = (b-a)/2^i\\)) and \\(j\\) is the
**extrapolation column**:

\\[
R[i, j] = \frac{4^j \, R[i, j-1] - R[i-1, j-1]}{4^j - 1}, \qquad 1 \leq j \leq i
\\]

The first column \\(R[i, 0]\\) is the composite trapezoidal estimate at step size
\\(h\_i\\). Column \\(j\\) contains estimates that are \\(O(h^{2(j+1)})\\) accurate. The
divisors \\(4^j - 1\\) take the values \\(3, 15, 63, 255, \ldots\\) as \\(j\\) increases.

For \\(n\\) columns the table has \\(n(n+1)/2\\) entries and the bottom-right entry
\\(R[n-1, n-1]\\) is the final estimate, with error \\(O(h^{2n})\\). For example, with
\\(n = 10\\) columns the error behaves as \\(O(h^{20})\\), which for smooth integrands
gives double-precision accuracy with far fewer evaluations than any fixed Newton-Cotes
rule of that order.

### Computation

**Midpoint reuse.** Each first-column entry \\(R[i, 0]\\) is obtained from the previous
one by adding only the \\(2^{i-1}\\) new midpoints — the \\(2^{i-1}+1\\) points already
evaluated at level \\(i-1\\) lie exactly on the level-\\(i\\) grid:

\\[
\begin{split}
R[i, 0] &= \frac{1}{2}\, R[i-1, 0] \\\\
&+ h\_i \sum\_{k=0}^{2^{i-1}-1} f\left(a + (2k+1)\,h\_i\right)
\end{split}
\\]

This means the total number of function evaluations across all \\(n\\) levels is
\\(2^{n-1}+1\\), not \\(O(n \cdot 2^{n-1})\\) as would be required if each level
recomputed all its points independently.

**Rolling buffer.** The Richardson extrapolation depends only on the previous row and
the current row of the table. The implementation therefore keeps only two rows of length
\\(n\\) in memory at any time, rather than the full \\(n \times n\\) table, keeping the
memory footprint at \\(O(n)\\).

**Summary of cost** for \\(n\\) columns:

| Quantity | Cost |
|---|---|
| Function evaluations | \\(2^{n-1}+1\\) |
| Memory | \\(O(n)\\) |
| Accuracy (smooth \\(f\\)) | \\(O(h^{2n})\\), \\(h = (b-a)/2^{n-1}\\) |

## Limitations

Romberg's method requires the integrand to be infinitely differentiable on \\([a, b]\\)
for the Euler-Maclaurin expansion to hold and the cancellation to be effective. It is not
appropriate for functions with singularities, discontinuities, or even kinks
(non-differentiable points) in the integration interval — at such points the error
expansion breaks down and the extrapolation accelerates to the wrong limit.

The number of columns \\(n\\) controls both the accuracy and the number of evaluations.
For double-precision arithmetic, \\(n \approx 10\\) is typically sufficient for smooth
integrands; using larger \\(n\\) beyond the point where \\(R[n-1, n-1]\\) has converged
to machine precision adds evaluations without improving the result. For functions with
near-singularities or rapid oscillation, prefer the adaptive Simpson method.
