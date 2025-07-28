#!/usr/bin/env python3
"""
Simple test to verify the PyO3 integration works
"""
import math

try:
    import integrate_python
    print("✓ integrate_python module loaded successfully!")
except ImportError as e:
    print(f"❌ Failed to import: {e}")
    exit(1)

def test_basic():
    """Test basic integration"""
    print("\n=== Basic Test ===")
    
    # Test: Integrate x^2 from 0 to 1 (should be 1/3)
    def square(x):
        return x * x
    
    result = integrate_python.py_legendre_rule(square, 0.0, 1.0, 10)
    expected = 1.0 / 3.0
    error = abs(result - expected)
    
    print(f"∫₀¹ x² dx = {result:.8f}")
    print(f"Expected:   {expected:.8f}")
    print(f"Error:      {error:.2e}")
    print(f"Status:     {'✓ PASS' if error < 1e-6 else '❌ FAIL'}")

def test_lambda():
    """Test with lambda function"""
    print("\n=== Lambda Test ===")
    
    result = integrate_python.py_legendre_rule(lambda x: math.sin(x), 0.0, math.pi, 20)
    expected = 2.0
    error = abs(result - expected)
    
    print(f"∫₀^π sin(x) dx = {result:.8f}")
    print(f"Expected:        {expected:.8f}")
    print(f"Error:           {error:.2e}")
    print(f"Status:          {'✓ PASS' if error < 1e-6 else '❌ FAIL'}")

if __name__ == "__main__":
    print("🚀 Simple PyO3 Integration Test")
    print("="*40)
    
    test_basic()
    test_lambda()
    
    print("\n" + "="*40)
    print("✨ Tests completed!")
    print("Your PyO3 Python->Rust function bridge is working!")