# Trapezoidal rule

The trapezoidal rule approximates the integral of a function \\(f(x)\\) on the closed and bounded interval \\(\[a, a+h\]\\)
of length \\(h > 0\\) by the (signed) area of the trapezoid formed by the line segments joining \\((a, 0)\\) to \\((a+h, 0)\\),
\\((a+h, 0)\\) to \\((a+h, f(a+h))\\),\\((a+h, f(a+h))\\) to \\((a, f(a))\\) and \\((a, f(a))\\) to \\((a, 0)\\).

The composite trapezoidal rule is used to approximate the integral of a function \\(f(x)\\) over a closed and bounded interval
\\(\[a, b\]\\) where \\(a < b\\), by decomposing the interval \\(\[a, b\]\\) into \\(n > 1\\) subintervals of equal length \\(h = \dfrac{b-a}{n}\\), then adding the results of applying the trapezoidal rule to each subinterval.

By abuse of language both the composite trapezoidal rule and the trapezoidal rule sometimes are referred to simply as
the trapezoidal rule.

Let \\(\int\_{a}^{b} f(x) dx\\) be the integral of \\(f(x)\\) over the closed and bounded interval \\(\[a,b\]\\), and
let \\(T_h(f)\\) be the result of applying the trapezoidal rule with \\(n\\) subintervals of length \\(h\\), i.e.

\\[
T_h(f) = h \left[ \frac{f(a)}{2} + f(a+h) + ··· + f(b-h) + \frac{f(b)}{2} \right]
\\]

The Euler-Maclaurin summation formula relates \\(\int\_{a}^{b} f(x) dx\\) and \\(T_h (f)\\)

\\[
\begin{align}
T_h(f) &= \int\_{a}^{b} f(x) dx + \frac{h^2}{12} \left[f'(b) - f'(a)\right] - \frac{h^4}{720} \left[f^{(3)}(b) - f^{(3)}(a) \right] \\\\
&+ ... + K h^{2p-2} \left[f^{(2p-3)}(b) - f^{(2p-3)(a)}\right] + O(h^{2p})
\end{align}
\\]

where \\(f^\prime\\), \\(f^{(3)}\\), and \\(f^{(2p-3)}\\) are the first, third and \\((p-3)^{th}\\) derivatives
of \\(f\\) and \\(K\\) is a constant.

The last term, \\(O(h^{2p})\\) is important. Given an infinitely differentiable function
in which the first \\((2p-3)\\) derivatives vanish at both endpoints of the interval of integration,
it is not true that \\(T_h (f) = \int\_{a}^{b} f(x) dx\\), but rather what the theorem says is that

\\[
\lim_{h \to 0} \left| \frac{T_h(f) - \int_{a}^{b} f(x)dx}{h^{2p}} \right| < M
\\]

where \\(M>0\\).

If \\(f\\) is at least twice differentiable on the interval \\(\[a,b\]\\), then applying the mean-value theorem to

\\[
\begin{align}
T_h(f) - \int\_{a}^{b} f(x) dx &= \frac{h^2}{12} \left[f^\prime(b) - f^\prime(a) \right] - \frac{h^4}{720} \left[ f^{(3)}(b) - f^{(3)}(a) \right] \\\\
&+ ... + K h^{2p-2} \left[f^{(2p-3)}(b) - f^{(2p-3)(a)} \right] + O(h^{2p})
\end{align}
\\]

yields the standard truncation error expression
\\[
R\_h(f) - \int\_{a}^{b} f(x) dx = -\frac{h^2}{12} (b - a) f^{\prime\prime}(c)
\\]

for some point \\(c\\) where \\(a ≤ c ≤ b\\).

A corollary of which is that if \\(f^{\prime\prime}(x) = 0\\) for all \\(x\\) in \\(\[a,b\]\\), i.e. if \\(f(x)\\) is linear,
then the trapezoidal rule is exact.

The Euler-Maclaurin summation formula also shows that usually \\(n\\) should be chosen large enough so that

\\[
h = \frac{b-a}{n} < 1
\\]

For example, if \\(h = 0.1\\) then
\\[
T\_{0.1}(f) = \int\_{a}^{b} f(x) dx - 0.00083 \left[f^\prime(b) - f^\prime(a) \right] + 0.00000014 \left[f^{(3)}(b) - f^{(3)}(a) \right] + ...
\\]
and if \\(h = 0.01\\) then
\\[
T\_{0.01}(f) = \int\_{a}^{b} f(x) dx - 0.0000083 \left[ f^\prime(b) - f^\prime(a) \right] + 0.000000000014 \left[f^{(3)}(b) - f^{(3)}(a)\right] + ...
\\]

while if \\(h=10\\) then
\\[
T\_{10}(f) = \int\_{a}^{b} f(x) dx - 8.3333 \left[ f^\prime(b) - f^\prime(a) \right] + 13.89 \left[ f^{(3)}(b) - f^{(3)}(a) \right] + ...
\\]

However, if the function \\(f(x)\\) is linear, then \\(n\\) may be chosen to be \\(1\\).
