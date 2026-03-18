# Adaptive Simpson

## Motivation

Fixed-step quadrature methods use the same subinterval size everywhere, which is inefficient when the integrand is smooth in most of the domain but has rapid variation or near-singularities in a small region. Adaptive quadrature automatically concentrates function evaluations where the integrand is hardest to approximate, achieving a target accuracy with fewer total evaluations.

## Example

```rust, editable
use integrate::adaptive_quadrature::adaptive_simpson_method;


let f = |x: f64| x.exp();

let a = 0.0;
let b = 1.0;

let tolerance = 10.0e-6;
let min_h = 10.0e-3;


let result = adaptive_simpson_method(f, a, b, min_h, tolerance);


match result{
    Ok(res)=>{println!("{}", res)}
    Err(err)=>{println!("{}", err)}
};
```

## Understanding adaptive Simpson method

If an integrand is poorly behaved in a small interval about a point,
then an attempt to integrate the function over an interval which contains
the poorly behaved interval either requires that small subintervals
are chosen for composite quadratures or the interval is decomposed into three intervals,
two on which the function is well-behaved and relatively large subintervals
can be chosen for the composite quadrature technique and one in which smaller subintervals need to be chosen.

Adaptive techniques are attempts to automatically detect and control the length of subintervals.

The technique for which the link to the listing is given below uses Simpson's rule
for integrating a function \\(f(x)\\) on a closed and bounded interval \\(\[a,b\]\\).

## Limitations

The adaptive Simpson method requires the integrand to be continuous on \\([a, b]\\). It may fail or produce unreliable results for functions with true singularities or jump discontinuities. The recursion depth is implicitly bounded by `min_h`; choosing `min_h` too large may miss fine structure, while too small a value can cause excessive recursion.
