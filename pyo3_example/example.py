#!/usr/bin/env python3
"""
Example usage of the integrate_python module showing how to pass
Python functions to Rust for fast numerical integration.
"""
import math
import time

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

# Import our Rust module (after building with maturin)
try:
    import integrate_python

    print("✓ integrate_python module loaded successfully!")
except ImportError:
    print("❌ integrate_python module not found. Please build it first with:")
    print("   cd pyo3_example && maturin develop")
    exit(1)


def test_basic_functions():
    """Test integration of basic mathematical functions."""
    print("\n=== Testing Basic Functions ===")

    # Test 1: Integrate x^2 from 0 to 1 (analytical result = 1/3)
    def square(x):
        return x * x

    result = integrate_python.py_legendre_rule(square, 0.0, 1.0, 100)
    analytical = 1.0 / 3.0
    error = abs(result - analytical)

    print(f"∫₀¹ x² dx:")
    print(f"  Numerical:  {result:.10f}")
    print(f"  Analytical: {analytical:.10f}")
    print(f"  Error:      {error:.2e}")
    print(f"  ✓ {'PASS' if error < 1e-10 else 'FAIL'}")

    # Test 2: Integrate sin(x) from 0 to π (analytical result = 2)
    def sin_func(x):
        return math.sin(x)

    result = integrate_python.py_legendre_rule(sin_func, 0.0, math.pi, 100)
    analytical = 2.0
    error = abs(result - analytical)

    print(f"\n∫₀^π sin(x) dx:")
    print(f"  Numerical:  {result:.10f}")
    print(f"  Analytical: {analytical:.10f}")
    print(f"  Error:      {error:.2e}")
    print(f"  ✓ {'PASS' if error < 1e-10 else 'FAIL'}")


def test_lambda_functions():
    """Test integration using Python lambda functions."""
    print("\n=== Testing Lambda Functions ===")

    # Test with lambda function
    result = integrate_python.py_legendre_rule(lambda x: x**3, 0.0, 2.0, 100)
    analytical = 2**4 / 4  # = 4
    error = abs(result - analytical)

    print(f"∫₀² x³ dx:")
    print(f"  Numerical:  {result:.10f}")
    print(f"  Analytical: {analytical:.10f}")
    print(f"  Error:      {error:.2e}")
    print(f"  ✓ {'PASS' if error < 1e-10 else 'FAIL'}")


def test_complex_functions():
    """Test integration of more complex functions."""
    print("\n=== Testing Complex Functions ===")

    # Test: Gaussian function (approximating ∫₋₃³ exp(-x²) dx)
    def gaussian(x):
        return math.exp(-x * x)

    result = integrate_python.py_legendre_rule(gaussian, -3.0, 3.0, 200)
    # Analytical result for ∫₋∞^∞ exp(-x²) dx = √π ≈ 1.7725,
    # but we're integrating from -3 to 3, which captures ~99.7% of the area
    expected_approx = math.sqrt(math.pi) * 0.997  # approximately

    print(f"∫₋₃³ exp(-x²) dx:")
    print(f"  Numerical:     {result:.8f}")
    print(f"  Expected ~:    {expected_approx:.8f}")
    print(f"  Ratio to √π:   {result/math.sqrt(math.pi):.6f}")


def test_performance():
    """Compare performance with different numbers of quadrature points."""
    print("\n=== Performance Test ===")

    def test_func(x):
        return math.sin(x) * math.exp(-x / 10)

    points = [10, 50, 100, 500, 1000]

    print("Points | Result      | Time (ms) | Error")
    print("-------|-------------|-----------|----------")

    for n in points:
        start = time.time()
        result = integrate_python.py_legendre_rule(test_func, 0.0, 10.0, n)
        elapsed = (time.time() - start) * 1000

        # Use the highest precision result as reference
        if n == 1000:
            reference = result
            error = 0.0
        else:
            error = abs(result - reference) if "reference" in locals() else 0.0

        print(f"{n:6d} | {result:11.8f} | {elapsed:9.3f} | {error:.2e}")


def test_numpy_integration():
    """Compare with NumPy if available."""
    if not HAS_NUMPY:
        print("\n=== NumPy Comparison (SKIPPED - NumPy not available) ===")
        return

    print("\n=== NumPy Comparison ===")
    from scipy import integrate as scipy_integrate

    def test_func(x):
        return x**2 * np.sin(x)

    # Our Rust implementation
    start = time.perf_counter()
    rust_result = integrate_python.py_legendre_rule(test_func, 0.0, math.pi, 8)
    rust_time = time.perf_counter() - start

    # SciPy quad for comparison
    start = time.perf_counter()
    scipy_result, scipy_error = scipy_integrate.quad(test_func, 0.0, math.pi)
    scipy_time = time.perf_counter() - start

    print(f"Function: x² sin(x) from 0 to π")
    print(f"Rust result:  {rust_result:.10f} ({rust_time*1000:.2f} ms)")
    print(f"SciPy result: {scipy_result:.10f} ({scipy_time*1000:.2f} ms)")
    print(f"Difference:   {abs(rust_result - scipy_result):.2e}")
    print(f"Speedup:      {scipy_time/rust_time:.1f}x")


def test_error_handling():
    """Test error handling with invalid inputs."""
    print("\n=== Error Handling ===")

    try:
        # This should raise an error (n=0)
        integrate_python.py_legendre_rule(lambda x: x, 0.0, 1.0, 0)
        print("❌ Should have raised an error for n=0")
    except ValueError as e:
        print(f"✓ Correctly caught error: {e}")

    # Test with a function that raises an exception
    def bad_func(x):
        if x > 0.5:
            raise ValueError("Test error")
        return x

    try:
        integrate_python.py_legendre_rule(bad_func, 0.0, 1.0, 10)
        print("❌ Should have raised an error for bad function")
    except:
        print("✓ Correctly handled function that raises exceptions")


if __name__ == "__main__":
    print("🚀 Testing PyO3 Integration Example")
    print("=" * 50)

    test_basic_functions()
    test_lambda_functions()
    test_complex_functions()
    test_performance()
    test_numpy_integration()
    # test_error_handling()

    print("\n" + "=" * 50)
    print("✨ All tests completed!")
    print(
        "\nThis demonstrates passing Python functions to Rust for fast numerical integration."
    )
    print("The PyFunction wrapper handles the Python->Rust function call bridge.")
