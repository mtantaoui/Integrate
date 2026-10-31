//! Integrand abstraction shared by the scalar and batch rules.
//!
//! Every fixed-node rule is written once, as a function generic over
//! [`Integrand`], and computes $\sum_i w_i f(x_i)$ through
//! [`Integrand::accumulate`]. The public scalar rules wrap the user's
//! `Fn(F1) -> F2` in [`Scalar`]; the `_batch` rules wrap the user's
//! `FnMut(&mut [F])` in [`Batch`].

use std::ops::Range;

use num::Float;

/// Maximum number of nodes passed to a batch integrand in a single call.
pub const CHUNK_SIZE: usize = 2048;

/// An integrand that can be evaluated on a range of indexed nodes.
pub trait Integrand<F: Float> {
    /// Adds $\sum_{i \in \text{range}} \text{weight}(i)_k \\, f(\text{node}(i))$
    /// to `sums[k]` for every weight set `k`, accumulating in `f64`.
    fn accumulate<const W: usize>(
        &mut self,
        range: Range<usize>,
        node: impl Fn(usize) -> F,
        weight: impl Fn(usize) -> [f64; W],
        sums: &mut [f64; W],
    );

    /// Convenience wrapper: $\sum_{i < \text{total}} w_i f(x_i)$ for each weight set.
    #[inline]
    fn weighted_sums<const W: usize>(
        &mut self,
        total: usize,
        node: impl Fn(usize) -> F,
        weight: impl Fn(usize) -> [f64; W],
    ) -> [f64; W] {
        let mut sums = [0.0_f64; W];
        self.accumulate(0..total, node, weight, &mut sums);
        sums
    }
}

#[inline]
fn add_weighted<const W: usize>(sums: &mut [f64; W], w: [f64; W], fx: f64) {
    for (s, wk) in sums.iter_mut().zip(w) {
        *s += wk * fx;
    }
}

/// Point-wise integrand `Fn(F1) -> F2`, evaluated one node at a time.
///
/// The value is converted straight to `f64`, so an `f32 -> f64` integrand
/// keeps its full output precision.
pub struct Scalar<Func>(pub Func);

impl<F1, F2, Func> Integrand<F1> for Scalar<Func>
where
    F1: Float,
    F2: Float,
    Func: Fn(F1) -> F2,
{
    #[inline]
    fn accumulate<const W: usize>(
        &mut self,
        range: Range<usize>,
        node: impl Fn(usize) -> F1,
        weight: impl Fn(usize) -> [f64; W],
        sums: &mut [f64; W],
    ) {
        for i in range {
            let fx = (self.0)(node(i))
                .to_f64()
                .expect("failed to convert f(x) to f64");
            add_weighted(sums, weight(i), fx);
        }
    }
}

/// Slice integrand `FnMut(&mut [F])` that overwrites each node with its value.
///
/// Nodes are evaluated in chunks of at most [`CHUNK_SIZE`] through a scratch
/// buffer that is allocated lazily and reused across calls.
pub struct Batch<F, Func> {
    func: Func,
    buf: Vec<F>,
}

impl<F, Func> Batch<F, Func> {
    pub fn new(func: Func) -> Self {
        Batch {
            func,
            buf: Vec::new(),
        }
    }
}

impl<F, Func> Integrand<F> for Batch<F, Func>
where
    F: Float,
    Func: FnMut(&mut [F]),
{
    fn accumulate<const W: usize>(
        &mut self,
        range: Range<usize>,
        node: impl Fn(usize) -> F,
        weight: impl Fn(usize) -> [f64; W],
        sums: &mut [f64; W],
    ) {
        let mut start = range.start;
        while start < range.end {
            let len = (range.end - start).min(CHUNK_SIZE);
            if self.buf.len() < len {
                self.buf.resize(len, F::zero());
            }
            let xs = &mut self.buf[..len];

            for (j, x) in xs.iter_mut().enumerate() {
                *x = node(start + j);
            }

            (self.func)(xs);

            for (j, &fx) in xs.iter().enumerate() {
                let fx = fx.to_f64().expect("failed to convert f(x) to f64");
                add_weighted(sums, weight(start + j), fx);
            }

            start += len;
        }
    }
}

/// `n * per + 1`, panicking with a clear message instead of overflowing.
pub fn node_count(n: usize, per: usize) -> usize {
    n.checked_mul(per)
        .and_then(|m| m.checked_add(1))
        .expect("number of intervals too large: node count overflows usize")
}
