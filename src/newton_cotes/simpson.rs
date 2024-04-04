//! Simpson's rule approximates the integral of a function $f(x)$ on the closed and bounded interval
//! $[a, a+h]$ of length $h > 0$ by the integral on $[a, a+h]$ of the quadratic passing through the
//! points $(a, f(a))$, $(a+h/2, f(a+h/2))$ and $(a+h, f(a+h))$. The composite Simpson's rule is used to
//! approximate the integral of a function $f(x)$ over a closed and bounded interval $[a, b]$ where $a < b$,
//! by decomposing the interval $[a, b]$ into $n > 1$ subintervals of equal length $h = \frac{b-a}{n}$,
//! then adding the results of applying the Simpson's rule to each subinterval. By abuse of
//! language both the composite Simpson's rule and Simpson's rule sometimes are referred to simply as
//! Simpson's rule. Let $int_{a}^{b} f(x) dx $ be the integral of $f(x)$ over the closed and bounded interval
//! $[a,b]$, and let $S_h(f)$ be the result of applying the Simpson's rule with $n$ subintervals of length h, i.e.
