# integrate-py

[![PyPI](https://img.shields.io/pypi/v/integrate-py)](https://pypi.org/project/integrate-py/)
[![Python](https://img.shields.io/pypi/pyversions/integrate-py)](https://pypi.org/project/integrate-py/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Python bindings for the [`integrate`](https://crates.io/crates/integrate) Rust crate.
Fast numerical integration in Rust, exposed as a Python package.

## Installation

```bash
pip install integrate-py
```

## Quick start

```python
import math
import integrate_py as ig

# ∫₀¹ x² dx = 1/3
ig.legendre_rule(lambda x: x**2, 0.0, 1.0, 100)

# ∫₀¹ eˣ dx = e − 1  (adaptive, tolerance 1e-6)
ig.adaptive_simpson(math.exp, 0.0, 1.0, 1e-3, 1e-6)

# ∫₀^π sin(x) dx = 2  (Romberg, 10 columns)
ig.romberg(math.sin, 0.0, math.pi, 10)
```

## Methods

| Function | Integral approximated |
|---|---|
| `rectangle_rule(f, a, b, n)` | ∫ₐᵇ f(x) dx — midpoint rectangle rule |
| `trapezoidal_rule(f, a, b, n)` | ∫ₐᵇ f(x) dx — trapezoidal rule |
| `simpson_rule(f, a, b, n)` | ∫ₐᵇ f(x) dx — Simpson 1/3 rule |
| `newton_rule(f, a, b, n)` | ∫ₐᵇ f(x) dx — Newton 3/8 rule |
| `legendre_rule(f, a, b, n)` | ∫ₐᵇ f(x) dx — Gauss-Legendre quadrature |
| `gauss_laguerre_rule(f, n)` | ∫₀^∞ f(x) e^{−x} dx |
| `gauss_hermite_rule(f, n)` | ∫_{−∞}^{+∞} f(x) e^{−x²} dx |
| `gauss_chebyshev1_rule(f, n)` | ∫_{−1}^{1} f(x) / √(1−x²) dx |
| `gauss_chebyshev2_rule(f, n)` | ∫_{−1}^{1} f(x) √(1−x²) dx |
| `adaptive_simpson(f, a, b, min_h, tol)` | ∫ₐᵇ f(x) dx — adaptive Simpson |
| `romberg(f, a, b, n_cols)` | ∫ₐᵇ f(x) dx — Romberg / Richardson extrapolation |

## Error handling

`adaptive_simpson` raises `ValueError` when it cannot meet the requested tolerance before
reaching `min_h`. All other methods raise `ValueError` for invalid arguments (e.g. `n=0`).

## Building from source

Requires Rust ≥ 1.63 and [maturin](https://www.maturin.rs/).

```bash
git clone https://github.com/mtantaoui/Integrate
cd Integrate/python
pip install maturin
maturin develop --release
```

## License

MIT — see [LICENSE](LICENSE).
