# Rectangle Rule

The rectangle rule approximates the integral of a function \\(f(x)\\) on the
closed and bounded interval \\([a, a+h]\\) of length \\(h > 0\\) by the (signed) area
of the rectangle with length \\(h\\) and height the value of the function \\(f(x)\\)
evaluated at the midpoint of the interval, \\(f(a + \frac{h}{2} )\\).

The composite rectangle rule is used to approximate the integral of a function
\\(f(x)\\) over a closed and bounded interval \\([a, b]\\) where \\(a < b\\), by decomposing
the interval \\([a, b]\\) into \\(n > 1\\) subintervals of equal length \\(h = \frac{b-a}{n}\\)
and adding the results of applying the rectangle rule to each subinterval.

By abuse of language both the composite rectangle rule and the rectangle rule sometimes
are referred to simply as the rectangle rule.

Let \\(\int\_{a}^{b} f(x) dx\\) be the integral of \\(f(x)\\) over the closed and bounded interval \\(\[a ,b \]\\),
and let \\(R_h(f)\\) be the result of applying the rectangle rule with \\(n\\) subintervals of length \\(h\\), i.e.

\\[
R_h(f)=h \left[ f(a+\frac{h}{2}) + f(a+\frac{3h}{2}) + ··· + f(b-\frac{h}{2}) \right]
\\]

An immediate consequence of the Euler-Maclaurin summation formula yields the following equation
relating \\(\int\_{a}^{b} f(x) dx\\) and \\(R_h(f)\\):

\begin{align}
R_h(f) & = \int\_{a}^{b} f(x) dx - \frac{h^2}{24} \left[ f^\prime (b) - f^\prime (a) \right] + \frac{7h^4}{5760} \left[ f^{(3)}(b) - f^{(3)}(a) \right] + \\\\ & ··· + K h^{2p-2} \left[ f^{(2p-3)}(b) - f^{(2p-3)}(a) \right] + O(h^{2p})\end{align}

where \\(f'\\), \\(f^{(3)}\\), and \\(f^{(2p-3)}\\) are the first, third and \\((2p-3)^{rd}\\) derivatives of \\(f\\) and \\(K\\) is a constant.

The last term, \\(O(h^{2p})\\) is important. Given an infinitely differentiable function
in which the first \\(2p-3\\) derivatives vanish at both endpoints of the interval of integration,
it is not true that \\(R*h(f) = \int*{a}^{b} f(x) dx\\), but rather what the theorem says is that

\\[
\lim_{h \to 0} \left| \dfrac{R_h(f) - \int_{a}^{b} f(x)dx}{h^{2p}} \right| < M
\\]

where \\(M>0\\).

If \\(f\\) is at least twice differentiable on the interval $\[a,b\]$, then applying the mean-value
theorem to

\begin{align} R_h(f) - \int\_{a}^{b} f(x) dx & = -\frac{h^2}{24} \left[ f^\prime (b) - f^\prime (a) \right] + \frac{7h^4}{5760} \left[ f^{(3)}(b) - f^{(3)}(a) \right] \\\\ & + ··· + K h^{2p-2} \left[ f^{(2p-3)}(b) - f^{(2p-3)}(a) \right] + O(h^{2p})\end{align}

yields the standard truncation error expression

\\[
R_h(f) - \int_{a}^{b} f(x) dx = -\frac{h^2}{24} (b - a) f^{\prime\prime}(c)
\\]

for some point \\(c\\) where \\(a ≤ c ≤ b\\).

A corollary of which is that if \\(f''(x) = 0\\) for all \\(x\\) in \\(\[a,b\]\\),
i.e. if \\(f(x)\\) is linear, then the rectangle rule is exact.

The Euler-Maclaurin summation formula also shows that usually \\(n\\) should be chosen large enough
so that \\(h = \frac{b-a}{n} < 1\\). For example, if \\(h = 0.1\\) then

\begin{split} R\_{0.1}(f) &= \int\_{a}^{b} f(x) dx - 0.00042 \left[ f'(b) - f'(a) \right] \\\\ & + 0.00000012 \left[f^{(3)}(b) - f^{(3)}(a) \right] + ... \end{split}

and if \\(h = 0.01\\) then

\begin{split} R\_{0.01}(f) &= \int\_{a}^{b} f(x) dx - 0.0000042 \left[ f^\prime(b) - f^\prime(a) \right] \\\\ &+ 0.000000000012 \left[ f^{(3)}(b) - f^{(3)}(a) \right] + ... \end{split}

while if \\(h=10\\) then

\begin{split} R\_{10}(f) &= \int\_{a}^{b} f(x) dx - 4.1667 \left[ f^\prime(b) - f^\prime(a)\right] \\\\ &+ 12.15 \left[ f^{(3)}(b) - f^{(3)}(a) \right] + ... \end{split}

However, if the function \\(f(x)\\) is linear, then $n$ may be chosen to be $1$.
