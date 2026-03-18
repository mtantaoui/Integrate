# Getting Started

## Installation

Add `integrate` to your `Cargo.toml`:

```toml
[dependencies]
integrate = "0.1"
```

## Your First Integral

Let's compute \\(\int_0^1 e^x \\, dx = e - 1 \approx 1.71828\\) using the trapezoidal rule:

```rust,editable
use integrate::newton_cotes::trapezoidal_rule;

let result = trapezoidal_rule(|x: f64| x.exp(), 0.0_f64, 1.0_f64, 1_000_usize);
println!("Result:   {:.6}", result);
println!("Exact:    {:.6}", std::f64::consts::E - 1.0);
```

## Using the Prelude

For convenience, you can import all methods at once with:

```rust
use integrate::prelude::*;
```

This exposes all 11 integration functions directly.

## Choosing a Method

| You want to integrate… | Recommended method |
|---|---|
| A smooth function on \\([a, b]\\) | [`legendre_rule`] or [`trapezoidal_rule`] |
| A function with variable smoothness | [`adaptive_simpson_method`] |
| A smooth function to high precision | [`romberg_method`] |
| \\(f(x)\\,e^{-x}\\) over \\([0, \infty)\\) | [`gauss_laguerre_rule`] |
| \\(f(x)\\,e^{-x^2}\\) over \\((-\infty, \infty)\\) | [`gauss_hermite_rule`] |
| \\(f(x)/\sqrt{1-x^2}\\) over \\([-1,1]\\) | [`gauss_chebyshev1_rule`] |
| \\(f(x)\sqrt{1-x^2}\\) over \\([-1,1]\\) | [`gauss_chebyshev2_rule`] |

## Generics: f32 and f64

All Newton-Cotes and Romberg methods are generic over float types. You can pass `f32` functions:

```rust,editable
use integrate::newton_cotes::simpson_rule;

let result = simpson_rule(|x: f32| x * x, 0.0_f32, 1.0_f32, 1_000_usize);
println!("{:.6}", result); // ≈ 0.333333
```
