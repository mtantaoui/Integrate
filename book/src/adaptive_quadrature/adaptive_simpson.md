# Adaptive Simpson

## Motivation

Fixed-step quadrature methods allocate function evaluations uniformly across the interval,
which is wasteful when the integrand is smooth almost everywhere but has rapid variation
or near-singularities in a small region. Adaptive quadrature automatically concentrates
evaluations where the integrand is hardest to approximate, achieving a target accuracy
with fewer total evaluations than any fixed-step rule of the same order.

The adaptive Simpson method combines the well-understood accuracy of Simpson's rule with
a local error indicator derived from Richardson extrapolation. Subintervals where the
integrand is smooth are accepted quickly; subintervals where it varies rapidly are
bisected and retried until the local error falls below a pro-rated fraction of the global
tolerance.

## Example

```rust,editable
use integrate::adaptive_quadrature::adaptive_simpson_method;

let f = |x: f64| x.exp();

let a = 0.0;
let b = 1.0;

let tolerance = 10.0e-6;
let min_h = 10.0e-3;

let result = adaptive_simpson_method(f, a, b, min_h, tolerance);

match result {
    Ok(res)  => println!("{}", res),
    Err(err) => println!("{}", err),
};
```

## Understanding the Adaptive Simpson method

### Error indicator

On a subinterval \\([l, r]\\) of length \\(H = r - l\\), two Simpson estimates are
computed from five evenly-spaced points \\(l,\\; l+H/4,\\; m,\\; r-H/4,\\; r\\)
(where \\(m = (l+r)/2\\)):

- **Single estimate** \\(S\_1\\) — standard Simpson using three points:

\\[
S\_1 = \frac{H}{6}\bigl[f(l) + 4f(m) + f(r)\bigr]
\\]

- **Composite estimate** \\(S\_2\\) — Simpson applied to the two half-panels
  \\([l, m]\\) and \\([m, r]\\):

\\[
S\_2 = \frac{H}{12}\bigl[f(l) + 4f(l+H/4) + 2f(m) + 4f(r-H/4) + f(r)\bigr]
\\]

By Richardson extrapolation, \\(S\_2\\) is approximately 16 times more accurate than
\\(S\_1\\) for smooth integrands (halving the step size reduces the \\(O(h^5)\\) per-panel
error by \\(2^4 = 16\\)). Consequently

\\[
S\_1 - \int\_{l}^{r} f\\,dx \approx \frac{16}{15}(S\_1 - S\_2)
\qquad \text{and} \qquad
\left|\int\_{l}^{r} f\\,dx - S\_2\right| \approx \frac{1}{15}|S\_1 - S\_2|
\\]

The difference \\(|S\_1 - S\_2|\\) is therefore a reliable local error indicator for
\\(S\_2\\), and the true error in the accepted estimate is roughly 15 times smaller than
the acceptance threshold — making the method conservative.

### Acceptance criterion

The global tolerance is distributed across subintervals in proportion to their length.
A subinterval \\([l, r]\\) is **accepted** when

\\[
|S\_1 - S\_2| < \varepsilon\_{\text{local}} = 2 \cdot \text{tol} \cdot \frac{r - l}{b - a}
\\]

where \\(\text{tol}\\) is the user-specified global tolerance and \\(b - a\\) is the full
integration interval. This proportional allocation ensures that the sum of all local
errors is bounded by the global tolerance.

### Algorithm

The implementation is **iterative**, using a linked stack of pending right-portions to
avoid recursion-depth limits that afflict the classical recursive formulation. Starting
from the full interval \\([a, b]\\):

1. Compute \\(S\_1\\) and \\(S\_2\\) on the current active subinterval.
2. If \\(|S\_1 - S\_2| < \varepsilon\_{\text{local}}\\): **accept** — add \\(S\_2\\) to the
   running total, pop the next pending subinterval from the stack, and reuse the already
   computed endpoint values \\(f(r)\\) and \\(f(m\_{\text{prev}})\\) as the new left
   endpoint and midpoint.
3. If \\(|S\_1 - S\_2| \geq \varepsilon\_{\text{local}}\\): **reject** — push the right
   half \\([m, r]\\) (with its cached function values) onto the stack, and set the active
   subinterval to the left half \\([l, m]\\).
4. Repeat from step 1. Terminate successfully when the right endpoint \\(b\\) is
   reached; terminate with an error if the active subinterval length falls below
   `min_h`.

Each accepted step reuses three of the five function values from the previous iteration
(the three points on the new coarser subinterval coincide with points already evaluated);
only two new evaluations (the quarter-points \\(l+H/4\\) and \\(r-H/4\\)) are needed per
step.

### Function evaluation budget

The number of evaluations depends on the integrand's behaviour:

- **Smooth integrands**: few bisections are needed; the algorithm accepts most
  subintervals on the first or second attempt, using roughly 2–4 evaluations per
  accepted panel beyond the initial three.
- **Rapidly varying integrands**: many bisections concentrate evaluations in the
  difficult region, while smooth regions are traversed cheaply.
- **Worst case**: if every step requires bisection down to `min_h`, the number of
  evaluations is proportional to \\((b-a)/\text{min\_h}\\).

## Limitations

The adaptive Simpson method requires the integrand to be continuous on \\([a, b]\\). It
may fail or produce unreliable results for functions with true singularities or jump
discontinuities inside the interval. The error indicator \\(|S\_1 - S\_2|\\) is only
reliable when the integrand is smooth on the current subinterval; it can be fooled by
integrand features that happen to cancel at the five sample points.

If a subinterval is bisected until its length reaches `min_h` without the tolerance
being satisfied, the method returns an error. This typically signals a singularity or
extremely rapid variation; try splitting the interval manually at the problem point,
reducing `min_h`, or increasing `tolerance`.
