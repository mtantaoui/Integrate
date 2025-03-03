# Adaptive Simpson

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
