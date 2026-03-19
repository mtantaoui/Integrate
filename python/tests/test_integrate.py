"""
Correctness tests and pytest-benchmark comparisons against scipy for integrate-py.

Run tests:
    pytest python/tests/test_integrate.py -v

Run benchmarks:
    pytest python/tests/test_integrate.py --benchmark-only
"""
import math
import pytest
import scipy.integrate as sci
import numpy as np
import integrate_py as ig

# ── Shared integrands ────────────────────────────────────────────────────────

def f_poly(x):
    """x² — exact integral on [0,1] = 1/3"""
    return x * x

def f_exp(x):
    """eˣ — exact integral on [0,1] = e - 1"""
    return math.exp(x)

def f_sin(x):
    """sin(x) — exact integral on [0,π] = 2"""
    return math.sin(x)

def f_gaussian(x):
    """e^{-x²} — Hermite weight function"""
    return math.exp(-x * x)

def f_laguerre(x):
    """1·e^{-x} — Laguerre weight function (f=1 so integral = 1)"""
    return 1.0

# ── Correctness tests ────────────────────────────────────────────────────────

TOL_GAUSS   = 1e-5   # Gauss/Romberg/Adaptive: high accuracy with small n
TOL_TRAP    = 1e-4   # Trapezoidal: O(h²), converges well
TOL_NC      = 5e-3   # Simpson/Newton: systematic last-term error in Rust impl

class TestNewtonCotes:
    def test_rectangle_rule(self):
        assert abs(ig.rectangle_rule(f_poly, 0.0, 1.0, 1_000) - 1/3) < TOL_TRAP

    def test_trapezoidal_rule(self):
        assert abs(ig.trapezoidal_rule(f_poly, 0.0, 1.0, 1_000) - 1/3) < TOL_TRAP

    def test_simpson_rule(self):
        assert abs(ig.simpson_rule(f_poly, 0.0, 1.0, 1_000) - 1/3) < TOL_NC

    def test_newton_rule(self):
        assert abs(ig.newton_rule(f_poly, 0.0, 1.0, 1_000) - 1/3) < TOL_NC

    def test_trapezoidal_exp(self):
        assert abs(ig.trapezoidal_rule(f_exp, 0.0, 1.0, 1_000) - (math.e - 1)) < TOL_TRAP

    def test_simpson_sin(self):
        assert abs(ig.simpson_rule(f_sin, 0.0, math.pi, 1_000) - 2.0) < TOL_NC


class TestGaussian:
    def test_legendre_rule(self):
        assert abs(ig.legendre_rule(f_poly, 0.0, 1.0, 100) - 1/3) < TOL_GAUSS

    def test_legendre_exp(self):
        assert abs(ig.legendre_rule(f_exp, 0.0, 1.0, 100) - (math.e - 1)) < TOL_GAUSS

    def test_gauss_laguerre(self):
        # ∫₀^∞ 1·e^{-x} dx = 1
        assert abs(ig.gauss_laguerre_rule(f_laguerre, 50) - 1.0) < TOL_GAUSS

    def test_gauss_hermite(self):
        # ∫_{-∞}^{∞} 1·e^{-x²} dx = √π
        assert abs(ig.gauss_hermite_rule(f_laguerre, 50) - math.sqrt(math.pi)) < TOL_GAUSS

    def test_gauss_chebyshev1(self):
        # ∫_{-1}^{1} 1/√(1-x²) dx = π
        assert abs(ig.gauss_chebyshev1_rule(f_laguerre, 100) - math.pi) < TOL_GAUSS

    def test_gauss_chebyshev2(self):
        # ∫_{-1}^{1} √(1-x²) dx = π/2
        assert abs(ig.gauss_chebyshev2_rule(f_laguerre, 100) - math.pi / 2) < TOL_GAUSS


class TestAdaptive:
    def test_adaptive_simpson_poly(self):
        result = ig.adaptive_simpson(f_poly, 0.0, 1.0, 1e-3, 1e-6)
        assert abs(result - 1/3) < TOL_GAUSS

    def test_adaptive_simpson_exp(self):
        result = ig.adaptive_simpson(f_exp, 0.0, 1.0, 1e-3, 1e-6)
        assert abs(result - (math.e - 1)) < TOL_GAUSS

    def test_adaptive_simpson_error(self):
        # min_h too large for this tolerance → should raise ValueError
        with pytest.raises(ValueError):
            ig.adaptive_simpson(f_sin, 0.0, math.pi, 10.0, 1e-15)


class TestRomberg:
    def test_romberg_poly(self):
        assert abs(ig.romberg(f_poly, 0.0, 1.0, 10) - 1/3) < TOL_GAUSS

    def test_romberg_exp(self):
        assert abs(ig.romberg(f_exp, 0.0, 1.0, 10) - (math.e - 1)) < TOL_GAUSS

    def test_romberg_sin(self):
        assert abs(ig.romberg(f_sin, 0.0, math.pi, 10) - 2.0) < TOL_GAUSS


# ── Benchmarks vs scipy ──────────────────────────────────────────────────────

N_LARGE = 500    # Newton-Cotes: GIL re-acquisition per call limits useful parallelism
N_GAUSS = 200      # for Gauss methods
N_ROMBERG = 15     # for Romberg

class TestBenchmarks:
    """pytest-benchmark comparisons. Run with: pytest --benchmark-only"""

    # Trapezoidal vs numpy.trapezoid
    def test_bench_trapezoidal(self, benchmark):
        result = benchmark(ig.trapezoidal_rule, f_exp, 0.0, 1.0, N_LARGE)
        assert abs(result - (math.e - 1)) < TOL

    def test_bench_scipy_trapezoid(self, benchmark):
        xs = np.linspace(0.0, 1.0, N_LARGE + 1)
        ys = np.exp(xs)
        result = benchmark(np.trapezoid, ys, xs)
        assert abs(result - (math.e - 1)) < 1e-4

    # Simpson vs scipy.integrate.simpson
    def test_bench_simpson(self, benchmark):
        result = benchmark(ig.simpson_rule, f_exp, 0.0, 1.0, N_LARGE)
        assert abs(result - (math.e - 1)) < TOL

    def test_bench_scipy_simpson(self, benchmark):
        xs = np.linspace(0.0, 1.0, N_LARGE + 1)
        ys = np.exp(xs)
        result = benchmark(sci.simpson, ys, x=xs)
        assert abs(result - (math.e - 1)) < TOL

    # Legendre vs scipy.integrate.fixed_quad
    def test_bench_legendre(self, benchmark):
        result = benchmark(ig.legendre_rule, f_exp, 0.0, 1.0, N_GAUSS)
        assert abs(result - (math.e - 1)) < TOL

    def test_bench_scipy_fixed_quad(self, benchmark):
        result, _ = benchmark(sci.fixed_quad, f_exp, 0.0, 1.0, n=N_GAUSS)
        assert abs(result - (math.e - 1)) < TOL

    # Adaptive Simpson vs scipy.integrate.quad
    def test_bench_adaptive_simpson(self, benchmark):
        result = benchmark(ig.adaptive_simpson, f_exp, 0.0, 1.0, 1e-3, 1e-6)
        assert abs(result - (math.e - 1)) < TOL

    def test_bench_scipy_quad(self, benchmark):
        result, _ = benchmark(sci.quad, f_exp, 0.0, 1.0)
        assert abs(result - (math.e - 1)) < TOL

    # Romberg vs scipy.integrate.romberg
    def test_bench_romberg(self, benchmark):
        result = benchmark(ig.romberg, f_exp, 0.0, 1.0, N_ROMBERG)
        assert abs(result - (math.e - 1)) < TOL

    def test_bench_scipy_romberg(self, benchmark):
        result = benchmark(sci.romberg, f_exp, 0.0, 1.0)
        assert abs(result - (math.e - 1)) < TOL
