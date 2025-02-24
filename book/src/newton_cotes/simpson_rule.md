# Simpson's rule

Simpson's rule approximates the integral of a function \\(f(x)\\) on the closed and bounded interval
\\(\[a, a+h\]\\) of length \\(h > 0\\) by the integral on \\(\[a, a+h\]\\) of the quadratic passing through the
points \\(\left(a, f(a)\right)\\), \\(\left(a+\dfrac{h}{2}, f\left(a+\dfrac{h}{2}\right)\right)\\) and \\(\left(a+h, f(a+h)\right)\\).

The composite Simpson's rule is used to approximate the integral of a function \\(f(x)\\) over a closed and bounded interval \\(\[a, b\]\\) where \\(a < b\\),
by decomposing the interval \\(\[a, b\]\\) into \\(n > 1\\) subintervals of equal length \\(h = \dfrac{b-a}{n}\\),
then adding the results of applying the Simpson's rule to each subinterval. By abuse of
language both the composite Simpson's rule and Simpson's rule sometimes are referred to simply as
Simpson's rule. Let \\(\int\_{a}^{b} f(x) dx\\) be the integral of \\(f(x)\\) over the closed and bounded interval
\\(\[a,b\]\\), and let \\(S_h(f)\\) be the result of applying the Simpson's rule with \\(n\\) subintervals of length h, i.e.

\\[
\begin{align}
S_h(f) &= \frac{h}{6} \left[ f(a) + 4f(a+\frac{h}{2}) + 2f(a+h) \right. \\\\
& + \left. ··· + 2f(b-h) + 4f(b-\frac{h}{2}) + f(b) \right]
\end{align}
\\]

An immediate consequence of the Euler-Maclaurin summation formula relates \\(\int\_{a}^{b} f(x) dx\\) and \\(S_h(f)\\):

\\[
\begin{align}
S_h(f) & = \int\_{a}^{b} f(x) dx + \frac{h^4}{2880} \left[ f^{(3)}(b) - f^{(3)}(a) \right] - \frac{h^6}{96768} \left[ f^{(5)}(b) - f^{(5)}(a) \right] \\\\
& + ··· + K h^{2p-2} \left[ f^{(2p-3)}(b) - f^{(2p-3)}(a) \right] + O(h^{2p})
\end{align}
\\]

where \\(f^{(3)}\\), \\(f^{(5)}\\), and \\(f^{(2p-3)}\\) are the third, fifth and \\((p-3)^{th}\\) derivatives of \\(f\\) and \\(K\\) is a constant.

The last term, \\(O(h^{2p})\\) is important. Given an infinitely differentiable function in which the first
\\(2p-3\\) derivatives vanish at both endpoints of the interval of integration, it is not true that

\\[
S_h(f) = \int_{a}^{b}f( x ) dx
\\]

but rather what the theorem says is that

\\[
\lim_{h \to 0} \left| \frac{S_h(f) - \int_{a}^{b} f(x)dx}{h^{2p}} \right| < M
\\]

where \\(M > 0\\).

If \\(f\\) is at least four times differentiable on the interval \\(\[a,b\]\\), then applying the mean-value theorem to

\\[
\begin{align}
S_h(f) - \int\_{a}^{b}f( x ) dx &= \frac{h^4}{2880} \left[ f^{(3)}(b) - f^{(3)}(a) \right] - \frac{h^6}{96768} \left[ f^{(5)}(b) - f^{(5)}(a) \right] \\\\
&+ ··· + K h^{(2p - 2)} \left[f^{(2p-3)}(b) - f^{(2p-3)}(a) \right] + O(h^{2p})
\end{align}
\\]

yields the standard truncation error expression

\\[
S_h(f) - \int_{a}^{b} f( x ) dx = \frac{h^4}{2880} (b-a) f^{(4)}(c)
\\]

for some point \\(c\\) where \\(a ≤ c ≤ b\\).

A corollary of which is that if \\(f^{(4)}(x) = 0\\) for all \\(x\\) in \\(\[a,b\]\\), i.e. if \\(f(x)\\) is a cubic,
then Simpson's rule is exact.

The Euler-Maclaurin summation formula also shows that usually $n$ should be chosen large enough so that \\(h = \frac{b-a}{n} < 1\\).

For example, if \\(h = 0.1\\) then
\\[
\begin{align}
S\_{0.1}(f) & = \int\_{a}^{b} f(x) dx + 3.5 · 10^{-8} \left[ f^{3}(b) - f^{3}(a) \right] \\\\
& - 1.033 · 10^{-11} \left[ f^{(5)}(b) - f^{(5)}(a) \right] + ···
\end{align}
\\]

and if \\(h = 0.01\\) then

\\[
\begin{split}
S\_{0.01}(f) & = \int\_{a}^{b} f(x) dx + 3.5 · 10^{-12} \left[ f^{3}(b) - f^{3}(a) \right] \\\\
& - 1.033 · 10^{-17} \left[ f^{(5)}(b) - f^{(5)}(a) \right] + ···
\end{split}
\\]

while if \\(h = 10\\) then
\\[
\begin{align}
S\_{0.01}(f) & = \int\_{a}^{b} f(x) dx + 3.47 \left[ f^{3}(b) - f^{3}(a) \right] \\\\
& - 10.33 \left[ f^{(5)}(b) - f^{(5)}(a) \right] + ···
\end{align}
\\]

However, if the function \\[f(x)\\] is a cubic, then \\[n\\] may be chosen to be \\[1\\].
