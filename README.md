# Integrate

![Integrate Logo](./book/src/images/Integrate.png)

<h3 style="text-align:center;">

[![Integrate crate](https://img.shields.io/crates/v/integrate.svg)](https://crates.io/crates/integrate)
[![Integrate documentation](https://img.shields.io/docsrs/integrate/latest)](https://docs.rs/integrate)
![minimum rustc 1.63](https://img.shields.io/badge/rustc-1.63+-red.svg)
[![build status](https://github.com/mtantaoui/Integrate/actions/workflows/main.yml/badge.svg)](https://github.com/mtantaoui/Integrate/actions)

</h3>

Integrate is a small, lightweight Rust library for performing numerical integration of real-valued functions. It is designed to integrate functions, providing a simple and efficient way to approximate definite integrals using various numerical methods.

Numerical integration is concerned with developing algorithms to
approximate the integral of a function $f(x)$. The most commonly used algorithms
are Newton-Cotes formulas, Romberg's method, Gaussian quadrature, and to
lesser extents Hermite's formulas and certain adaptive techniques.

## Features

Integrate supports a variety of numerical integration techniques:

- Newton-Cotes methods:

  - Rectangle Rule.
  - Trapezoidal Rule.
  - Simpson's Rule.
  - Newton's 3/8 Rule.

- Gauss quadrature methods:

  - Gauss-Legendre.
  - Gauss-Laguerre.
  - Gauss-Hermite.
  - Gauss-Chebyshev First Kind.
  - Gauss-Chebyshev Second Kind.

- Adaptive Methods:

  - Adaptive Simpson's method

- Romberg’s method.

## Caveats

All of the numerical integration techniques listed above assume that
the integrand is continuous on the interval of integration. The error
estimates generally require that the integrands are differentiable of
whatever order is required so that the formula for the error estimate
makes sense. Below is a check list which one should verify before using
any of the numerical algorithms.

### Integrand check list

- Are there any singularities of the integrand in the interval of integration?

  - If there are, then can they be removed?

  - A function which has jump discontinuities can be integrated by splitting
    the interval of integration into subintervals on which the function is continuous,
    and then numerically integrate over each subinterval and add the results.

  - Using integration by parts certain types of singular integrals can be
    expressed as the sum of a function and a non-singular integral.

  - In other cases, an approximation of the integral in a small neighborhood of
    the singularity can be obtained by hand and the interval of integration can be
    split into three intervals in which numerical integration can be used on
    two of them.

- Does the function oscillate over the region of integration? If it does,
  then make sure that the step size is chosen to be smaller than the wave length
  of the function. The interval of integration can also be split into subintervals
  in which each subinterval is half a wave length and the algorithm is applied
  to each subinterval.

### Installation

To use the `Integrate` crate in your Rust project, add the following line to your `Cargo.toml`:

```toml
[dependencies]
integrate = "0.1.7"
```

### Contribution

Feel free to submit issues or pull requests for bug fixes, new features, or other improvements. Contributions are welcome!

### License

This project is licensed under the MIT License – see the LICENSE file for details. Opening a pull request is assumed to signal agreement with these licensing terms.
