"""
integrate-py: Fast numerical integration via Rust.

All 11 methods are parallelized with Rayon under the hood.
Import with::

    from integrate_py import rectangle_rule, legendre_rule, ...

or use the module directly::

    import integrate_py as ig
    ig.trapezoidal_rule(lambda x: x**2, 0.0, 1.0, 1000)
"""

from ._integrate_py import (
    rectangle_rule,
    trapezoidal_rule,
    simpson_rule,
    newton_rule,
    legendre_rule,
    gauss_laguerre_rule,
    gauss_hermite_rule,
    gauss_chebyshev1_rule,
    gauss_chebyshev2_rule,
    adaptive_simpson,
    romberg,
)

__all__ = [
    "rectangle_rule",
    "trapezoidal_rule",
    "simpson_rule",
    "newton_rule",
    "legendre_rule",
    "gauss_laguerre_rule",
    "gauss_hermite_rule",
    "gauss_chebyshev1_rule",
    "gauss_chebyshev2_rule",
    "adaptive_simpson",
    "romberg",
]
