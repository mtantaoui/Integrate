# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed
- Remove Rayon parallelism dependency (sequential iterators perform better for typical workloads)
- Rewrite Romberg integration: iterative rolling-buffer algorithm, O(2^n) function evaluations
- Fix Newton's 3/8 rule: corrected factor and last-term evaluation bugs
- Fix Simpson's rule: last term now correctly evaluates f(b)
- Raise MSRV from 1.63 to 1.89 (required by the `simdly` dev-dependency used in tests and benchmarks)
- Every fixed-node rule is now implemented once and shared by its scalar and `_batch` variants; scalar results may differ from 0.3.x in the last bits (summation order)
- Gauss-Laguerre, Gauss-Hermite and Gauss-Chebyshev rules now accumulate in `f64` for every float type (more accurate for `f32`)
- `romberg_method` panics if `n_columns` exceeds `MAX_ROMBERG_COLUMNS` (32) instead of overflowing

### Added
- `batch` module: `_batch` variants of every fixed-node rule (rectangle, trapezoidal, Simpson, Newton 3/8, Romberg, Gauss-Legendre/Laguerre/Hermite/Chebyshev, Gauss-Kronrod) take an integrand that evaluates a whole slice of nodes in place (e.g. with SIMD libraries such as `simdly`)
- `romberg::MAX_ROMBERG_COLUMNS` constant
- Full rustdoc for all 11 public integration functions
- mdBook documentation with getting started guide and Python API chapter
- Python bindings via PyO3 + Maturin (`integrate_py` package)
- CI: multi-platform test matrix, lint, docs, and Pages deploy

### Fixed
- Simpson, Newton 3/8 and trapezoidal rules panic with a clear message instead of overflowing when the node count exceeds `usize::MAX`
- LaTeX rendering in KaTeX and MathJax (escaped chars, multi-line formulas)
- Wrong sign in Gauss-Hermite weight function documentation
- Wrong integration limits in Gauss-Chebyshev documentation
