# Python API

`integrate-py` exposes all 11 integration methods as a Python package, built with
[PyO3](https://pyo3.rs) and [Maturin](https://www.maturin.rs). The Rust core is called
directly from Python, with the GIL released during integration so other Python threads
can run concurrently.

## Installation

```bash
pip install integrate-py
```

Or build from source (requires Rust ≥ 1.63):

```bash
git clone https://github.com/mtantaoui/Integrate
cd Integrate/python
pip install maturin
maturin develop --release
```

## Quick start

```python
import math
import integrate_py as ig

# ∫₀¹ x² dx = 1/3
result = ig.legendre_rule(lambda x: x**2, 0.0, 1.0, 100)
print(f"{result:.10f}")   # 0.3333333333

# ∫₀¹ eˣ dx = e − 1
result = ig.adaptive_simpson(math.exp, 0.0, 1.0, min_h=1e-3, tol=1e-6)
print(f"{result:.10f}")   # 1.7182818284

# ∫₀^π sin(x) dx = 2
result = ig.romberg(math.sin, 0.0, math.pi, n_cols=10)
print(f"{result:.10f}")   # 2.0000000000
```

## Available methods

- [Newton-Cotes](python_api/newton_cotes.md)
- [Gaussian Quadrature](python_api/gaussian.md)
- [Adaptive Simpson](python_api/adaptive.md)
- [Romberg](python_api/romberg.md)

## Choosing a method

| You want to integrate… | Recommended |
|---|---|
| Smooth \\(f\\) on \\([a, b]\\) | `legendre_rule` or `romberg` |
| Highly smooth \\(f\\), maximum precision | `romberg` |
| \\(f\\) with variable smoothness | `adaptive_simpson` |
| Simple estimate, any \\(f\\) | `trapezoidal_rule` |
| \\(\int_0^\infty f(x)\\,e^{-x}\\,dx\\) | `gauss_laguerre_rule` |
| \\(\int_{-\infty}^{\infty} f(x)\\,e^{-x^2}\\,dx\\) | `gauss_hermite_rule` |
| \\(\int_{-1}^{1} f(x)/\sqrt{1-x^2}\\,dx\\) | `gauss_chebyshev1_rule` |
| \\(\int_{-1}^{1} f(x)\sqrt{1-x^2}\\,dx\\) | `gauss_chebyshev2_rule` |

## A note on performance

The Rust core runs without Python overhead between function evaluations. For best performance:

- **Prefer methods with small \\(n\\)** (`legendre_rule`, `romberg`, `adaptive_simpson`) — they
  need far fewer function evaluations than Newton-Cotes for the same accuracy.
- **Use NumPy-vectorized integrands** when very large \\(n\\) is needed; pre-evaluate into an
  array and integrate with `numpy.trapezoid` / `scipy.integrate.simpson`.
- **The GIL is released** during integration, so other Python threads are not blocked
  while the Rust core computes.

## Benchmarks vs scipy

The table below compares `integrate-py` against the equivalent
[`scipy.integrate`](https://docs.scipy.org/doc/scipy/reference/integrate.html) methods on
\\(\int_0^1 e^x\\,dx = e - 1\\), using `pytest-benchmark`. Results will vary by machine.

Run the benchmarks yourself:

```bash
cd Integrate
pip install pytest pytest-benchmark scipy numpy maturin
maturin develop --release --manifest-path python/Cargo.toml
pytest python/tests/test_integrate.py --benchmark-only -v
```

| Method | integrate-py | scipy equivalent |
|---|---|---|
| `trapezoidal_rule(n=500)` | `pytest-benchmark` | `numpy.trapezoid` |
| `simpson_rule(n=500)` | `pytest-benchmark` | `scipy.integrate.simpson` |
| `legendre_rule(n=200)` | `pytest-benchmark` | `scipy.integrate.fixed_quad` |
| `adaptive_simpson` | `pytest-benchmark` | `scipy.integrate.quad` |
| `romberg(n_cols=15)` | `pytest-benchmark` | `scipy.integrate.romberg` |
