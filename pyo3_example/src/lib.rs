use integrate::prelude::*;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// A wrapper struct to hold Python callable objects
/// This allows us to pass Python functions to Rust integration methods
struct PyFunction {
    func: PyObject,
}

impl PyFunction {
    fn new(func: PyObject) -> Self {
        Self { func }
    }

    /// Call the Python function with a single f64 argument
    fn call(&self, x: f64) -> PyResult<f64> {
        Python::with_gil(|py| {
            let result = self.func.call1(py, (x,))?;
            result.extract::<f64>(py)
        })
    }
}

/// Wrapper for Gauss-Legendre quadrature
///
/// Integrates a Python function using Gauss-Legendre quadrature
///
/// Args:
///     func: Python callable that takes a float and returns a float
///     a: Lower integration limit
///     b: Upper integration limit  
///     n: Number of quadrature points
///
/// Returns:
///     Numerical integral approximation
///
/// Example:
///     >>> import math
///     >>> def f(x):
///     ...     return x * x
///     >>> result = py_legendre_rule(f, 0.0, 1.0, 100)
///     >>> print(f"Integral of x^2 from 0 to 1: {result}")
#[pyfunction]
fn py_legendre_rule(py: Python, func: PyObject, a: f64, b: f64, n: usize) -> PyResult<f64> {
    if n == 0 {
        return Err(PyValueError::new_err(
            "Number of points must be greater than 0",
        ));
    }

    let py_func = PyFunction::new(func);

    // Create a closure that captures the Python function
    let rust_func = |x: f64| -> f64 {
        match py_func.call(x) {
            Ok(result) => result,
            Err(err) => {
                // Convert Python error to panic for now
                // In production, you might want better error handling
                panic!("Error calling Python function: {:?}", err);
            }
        }
    };

    py.allow_threads(|| Ok(legendre_rule(rust_func, a, b, n)))
}

/// Python module definition
#[pymodule]
fn integrate_python(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_legendre_rule, m)?)?;

    // Add module documentation
    m.add(
        "__doc__",
        "Fast Gauss-Legendre numerical integration with Python bindings",
    )?;

    Ok(())
}
