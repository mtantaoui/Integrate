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
- Raise MSRV from 1.63 to 1.71

### Added
- Full rustdoc for all 11 public integration functions
- mdBook documentation with getting started guide and Python API chapter
- Python bindings via PyO3 + Maturin (`integrate_py` package)
- CI: multi-platform test matrix, lint, docs, and Pages deploy

### Fixed
- LaTeX rendering in KaTeX and MathJax (escaped chars, multi-line formulas)
- Wrong sign in Gauss-Hermite weight function documentation
- Wrong integration limits in Gauss-Chebyshev documentation
