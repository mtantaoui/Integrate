use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::cell::RefCell;

// ── PyFunction helper ─────────────────────────────────────────────────────────

/// Wraps a Python callable. `PyObject` is `Send + Sync`, so `PyFunction` is too.
struct PyFunction(PyObject);

impl PyFunction {
    fn call(&self, x: f64) -> f64 {
        Python::with_gil(|py| {
            self.0
                .call1(py, (x,))
                .and_then(|r| r.extract::<f64>(py))
                .unwrap_or_else(|e| panic!("Python function call failed: {e}"))
        })
    }
}

// ── Thread-local for adaptive_simpson (Copy bound, single-threaded) ──────────

// `adaptive_simpson_method` requires `Func: Fn(F) -> F + Sync + Copy`.
// It is single-threaded/recursive, so a thread-local + fn pointer works.
thread_local! {
    static TL_FUNC: RefCell<Option<PyObject>> = const { RefCell::new(None) };
}

fn tl_call(x: f64) -> f64 {
    TL_FUNC.with(|slot| {
        Python::with_gil(|py| {
            slot.borrow()
                .as_ref()
                .expect("TL_FUNC not set")
                .call1(py, (x,))
                .and_then(|r| r.extract::<f64>(py))
                .unwrap_or_else(|e| panic!("Python function call failed: {e}"))
        })
    })
}

// ── Raw-pointer wrapper for romberg (Copy + Send + Sync bounds) ──────────────

// `romberg_method` requires `Func: Fn(F1) -> F2 + Send + Sync + Copy`.
// We wrap a raw pointer to the `PyObject` that lives on the caller's stack.
// Safety: the `PyObject` is valid for the entire duration of `romberg_method`
// (it's a local variable in the calling function), and each call re-acquires
// the GIL before touching Python objects.
#[derive(Clone, Copy)]
struct RawPyFunc(*const PyObject);

// SAFETY: `PyObject` is `Send + Sync`; we only call it while holding the GIL.
unsafe impl Send for RawPyFunc {}
unsafe impl Sync for RawPyFunc {}

impl RawPyFunc {
    fn call(self, x: f64) -> f64 {
        Python::with_gil(|py| {
            // SAFETY: pointer is valid for the duration of the enclosing fn.
            unsafe { &*self.0 }
                .call1(py, (x,))
                .and_then(|r| r.extract::<f64>(py))
                .unwrap_or_else(|e| panic!("Python function call failed: {e}"))
        })
    }
}

// ── Newton-Cotes ──────────────────────────────────────────────────────────────

#[pyfunction]
#[pyo3(text_signature = "(func, a, b, n, /)")]
fn rectangle_rule(py: Python, func: PyObject, a: f64, b: f64, n: usize) -> PyResult<f64> {
    if n == 0 {
        return Err(PyValueError::new_err("n must be > 0"));
    }
    let f = PyFunction(func);
    Ok(py.allow_threads(move || {
        integrate::newton_cotes::rectangle_rule(|x: f64| f.call(x), a, b, n)
    }))
}

#[pyfunction]
#[pyo3(text_signature = "(func, a, b, n, /)")]
fn trapezoidal_rule(py: Python, func: PyObject, a: f64, b: f64, n: usize) -> PyResult<f64> {
    if n == 0 {
        return Err(PyValueError::new_err("n must be > 0"));
    }
    let f = PyFunction(func);
    Ok(py.allow_threads(move || {
        integrate::newton_cotes::trapezoidal_rule(|x: f64| f.call(x), a, b, n)
    }))
}

#[pyfunction]
#[pyo3(text_signature = "(func, a, b, n, /)")]
fn simpson_rule(py: Python, func: PyObject, a: f64, b: f64, n: usize) -> PyResult<f64> {
    if n == 0 {
        return Err(PyValueError::new_err("n must be > 0"));
    }
    let f = PyFunction(func);
    Ok(py.allow_threads(move || {
        integrate::newton_cotes::simpson_rule(|x: f64| f.call(x), a, b, n)
    }))
}

#[pyfunction]
#[pyo3(text_signature = "(func, a, b, n, /)")]
fn newton_rule(py: Python, func: PyObject, a: f64, b: f64, n: usize) -> PyResult<f64> {
    if n == 0 {
        return Err(PyValueError::new_err("n must be > 0"));
    }
    let f = PyFunction(func);
    Ok(py.allow_threads(move || {
        integrate::newton_cotes::newton_rule(|x: f64| f.call(x), a, b, n)
    }))
}

// ── Gaussian quadrature ───────────────────────────────────────────────────────

#[pyfunction]
#[pyo3(text_signature = "(func, a, b, n, /)")]
fn legendre_rule(py: Python, func: PyObject, a: f64, b: f64, n: usize) -> PyResult<f64> {
    if n == 0 {
        return Err(PyValueError::new_err("n must be > 0"));
    }
    let f = PyFunction(func);
    Ok(py.allow_threads(move || {
        integrate::gauss_quadrature::legendre_rule(|x: f64| f.call(x), a, b, n)
    }))
}

#[pyfunction]
#[pyo3(text_signature = "(func, n, /)")]
fn gauss_laguerre_rule(py: Python, func: PyObject, n: usize) -> PyResult<f64> {
    if n == 0 {
        return Err(PyValueError::new_err("n must be > 0"));
    }
    let f = PyFunction(func);
    Ok(py.allow_threads(move || {
        integrate::gauss_quadrature::gauss_laguerre_rule::<_, f64>(|x: f64| f.call(x), n)
    }))
}

#[pyfunction]
#[pyo3(text_signature = "(func, n, /)")]
fn gauss_hermite_rule(py: Python, func: PyObject, n: usize) -> PyResult<f64> {
    if n == 0 {
        return Err(PyValueError::new_err("n must be > 0"));
    }
    let f = PyFunction(func);
    Ok(py.allow_threads(move || {
        integrate::gauss_quadrature::gauss_hermite_rule::<_, f64>(|x: f64| f.call(x), n)
    }))
}

#[pyfunction]
#[pyo3(text_signature = "(func, n, /)")]
fn gauss_chebyshev1_rule(py: Python, func: PyObject, n: usize) -> PyResult<f64> {
    if n == 0 {
        return Err(PyValueError::new_err("n must be > 0"));
    }
    let f = PyFunction(func);
    Ok(py.allow_threads(move || {
        integrate::gauss_quadrature::gauss_first_kind_chebyshev_rule::<_, f64>(
            |x: f64| f.call(x),
            n,
        )
    }))
}

#[pyfunction]
#[pyo3(text_signature = "(func, n, /)")]
fn gauss_chebyshev2_rule(py: Python, func: PyObject, n: usize) -> PyResult<f64> {
    if n == 0 {
        return Err(PyValueError::new_err("n must be > 0"));
    }
    let f = PyFunction(func);
    Ok(py.allow_threads(move || {
        integrate::gauss_quadrature::gauss_second_kind_chebyshev_rule::<_, f64>(
            |x: f64| f.call(x),
            n,
        )
    }))
}

// ── Adaptive quadrature ───────────────────────────────────────────────────────

// `adaptive_simpson_method` requires `Func: Fn(F) -> F + Sync + Copy`.
// We use a thread-local + fn pointer to satisfy the `Copy` bound.
#[pyfunction]
#[pyo3(text_signature = "(func, a, b, min_h, tol, /)")]
fn adaptive_simpson(
    _py: Python,
    func: PyObject,
    a: f64,
    b: f64,
    min_h: f64,
    tol: f64,
) -> PyResult<f64> {
    TL_FUNC.with(|slot| *slot.borrow_mut() = Some(func));
    let result = integrate::adaptive_quadrature::adaptive_simpson_method(
        tl_call as fn(f64) -> f64,
        a,
        b,
        min_h,
        tol,
    )
    .map_err(|e| PyValueError::new_err(e.to_string()));
    TL_FUNC.with(|slot| *slot.borrow_mut() = None);
    result
}

// ── Romberg ───────────────────────────────────────────────────────────────────

// `romberg_method` requires `Func: Fn(F1) -> F2 + Sync + Send + Copy`.
// Uses RawPyFunc — a Copy + Send + Sync wrapper around a raw pointer to `func`.
// `func` lives on the stack here, so the pointer is valid for the call duration.
#[pyfunction]
#[pyo3(text_signature = "(func, a, b, n_cols, /)")]
fn romberg(py: Python, func: PyObject, a: f64, b: f64, n_cols: usize) -> PyResult<f64> {
    if n_cols == 0 {
        return Err(PyValueError::new_err("n_cols must be > 0"));
    }
    let fp = RawPyFunc(&func as *const PyObject);
    let result = py.allow_threads(|| {
        integrate::romberg::romberg_method(move |x: f64| fp.call(x), a, b, n_cols)
    });
    Ok(result)
}

// ── Module ────────────────────────────────────────────────────────────────────

#[pymodule]
fn _integrate_py(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(rectangle_rule, m)?)?;
    m.add_function(wrap_pyfunction!(trapezoidal_rule, m)?)?;
    m.add_function(wrap_pyfunction!(simpson_rule, m)?)?;
    m.add_function(wrap_pyfunction!(newton_rule, m)?)?;
    m.add_function(wrap_pyfunction!(legendre_rule, m)?)?;
    m.add_function(wrap_pyfunction!(gauss_laguerre_rule, m)?)?;
    m.add_function(wrap_pyfunction!(gauss_hermite_rule, m)?)?;
    m.add_function(wrap_pyfunction!(gauss_chebyshev1_rule, m)?)?;
    m.add_function(wrap_pyfunction!(gauss_chebyshev2_rule, m)?)?;
    m.add_function(wrap_pyfunction!(adaptive_simpson, m)?)?;
    m.add_function(wrap_pyfunction!(romberg, m)?)?;
    Ok(())
}
