# PyO3 Integration Example

This example demonstrates how to use PyO3 to pass Python functions to your Rust numerical integration library. It focuses on the Gauss-Legendre quadrature method as a practical example.

## Key Concepts

### The Challenge
Passing Python functions to Rust is tricky because:
- Python functions are objects in the Python runtime
- Rust needs a callable that matches its function signature
- We need to bridge between Python's dynamic typing and Rust's static typing

### The Solution
The example uses a `PyFunction` wrapper struct that:
1. Holds a reference to the Python callable (`PyObject`)
2. Provides a `call()` method that invokes the Python function from Rust
3. Handles the Python GIL (Global Interpreter Lock) correctly
4. Converts between Python and Rust types

## Files

- `Cargo.toml` - Rust dependencies including PyO3 and your integrate library
- `pyproject.toml` - Python build configuration using Maturin
- `src/lib.rs` - Rust code with PyO3 bindings
- `example.py` - Python examples showing usage
- `README.md` - This documentation

## How It Works

### 1. PyFunction Wrapper
```rust
#[derive(Clone)]
struct PyFunction {
    func: PyObject,  // Holds the Python function
}

impl PyFunction {
    fn call(&self, x: f64) -> PyResult<f64> {
        Python::with_gil(|py| {
            let result = self.func.call1(py, (x,))?;
            result.extract::<f64>(py)
        })
    }
}
```

### 2. Rust Integration Function
```rust
#[pyfunction]
fn py_legendre_rule(func: PyObject, a: f64, b: f64, n: usize) -> PyResult<f64> {
    let py_func = PyFunction::new(func);
    
    let rust_func = |x: f64| -> f64 {
        match py_func.call(x) {
            Ok(result) => result,
            Err(_) => panic!("Error calling Python function"),
        }
    };
    
    Ok(legendre_rule(rust_func, a, b, n))
}
```

### 3. Python Usage
```python
import integrate_python

def my_function(x):
    return x * x

result = integrate_python.py_legendre_rule(my_function, 0.0, 1.0, 100)
```

## Build Instructions

### Prerequisites
```bash
# Install Rust
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Install Python development headers (Ubuntu/Debian)
sudo apt install python3-dev

# Install Maturin
pip install maturin
```

### Building
```bash
cd pyo3_example

# Development build (creates wheel and installs in current Python environment)
maturin develop

# Production build (creates wheel file)
maturin build --release
```

### Running the Example
```bash
python example.py
```

## Expected Output
```
🚀 Testing PyO3 Integration Example
==================================================

=== Testing Basic Functions ===
∫₀¹ x² dx:
  Numerical:  0.3333333333
  Analytical: 0.3333333333
  Error:      1.39e-11
  ✓ PASS

∫₀^π sin(x) dx:
  Numerical:  2.0000000000
  Analytical: 2.0000000000
  Error:      4.44e-16
  ✓ PASS

... (more tests)
```

## Performance Notes

- The PyO3 bridge has some overhead for each function call
- For compute-intensive integrands, this overhead is negligible
- For simple functions, pure Python implementations might be faster
- The Gauss-Legendre method typically needs fewer function evaluations than Newton-Cotes methods

## Common Issues & Solutions

### Issue: "ImportError: No module named 'integrate_python'"
**Solution:** Run `maturin develop` to build and install the module.

### Issue: "Python.h not found"
**Solution:** Install Python development headers:
- Ubuntu/Debian: `sudo apt install python3-dev`
- RedHat/CentOS: `sudo yum install python3-devel`
- macOS: Usually included with Python installation

### Issue: Function panics with "Error calling Python function"
**Solution:** Check that your Python function:
- Takes a single float argument
- Returns a number (int or float)
- Doesn't raise exceptions for the integration domain

## Advanced Usage

### Error Handling
```python
try:
    result = integrate_python.py_legendre_rule(bad_function, 0.0, 1.0, 100)
except ValueError as e:
    print(f"Integration failed: {e}")
```

### NumPy Arrays
```python
import numpy as np

def vectorized_func(x):
    # Note: x is a scalar here, not an array
    return np.sin(x) * np.exp(-x)

result = integrate_python.py_legendre_rule(vectorized_func, 0.0, 10.0, 100)
```

### Lambda Functions
```python
result = integrate_python.py_legendre_rule(lambda x: x**3 + 2*x, -1.0, 1.0, 50)
```

## Extension Ideas

1. **Add more integration methods** - Extend to other methods from your library
2. **Vectorized functions** - Support functions that take numpy arrays
3. **Complex numbers** - Support complex-valued functions
4. **Parallel integration** - Integrate multiple functions simultaneously
5. **Adaptive methods** - Expose adaptive Simpson and other error-controlled methods

## Technical Details

- Uses PyO3 0.22 for Python-Rust bindings
- Requires Python 3.8+
- Built with Maturin for easy Python packaging
- Handles Python GIL automatically
- Provides proper error propagation from Python to Rust

This example provides a solid foundation for exposing your high-performance Rust numerical integration library to Python users!