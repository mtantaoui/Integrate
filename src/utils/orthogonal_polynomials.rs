use std::{fmt::Debug, ops::AddAssign};

use num::Float;

/// A family of polynomials that are mutually orthogonal with respect to a
/// weight function on a given interval.
///
/// Implementations of this trait are used internally to compute the nodes
/// (zeros of the polynomial) and weights required by Gaussian quadrature rules.
/// The nodes are derived from the eigenvalues of the associated Jacobi
/// tridiagonal matrix via the Golub-Welsch algorithm.
pub trait OrthogonalPolynomial<F: Float + Debug + AddAssign> {
    /// Creates a polynomial of the given `degree`.
    fn new(degree: usize) -> Self;

    /// Evaluates the polynomial at `x` using the three-term recurrence relation.
    fn eval(&self, x: F) -> F;

    /// Returns the zeros of the polynomial.
    ///
    /// The zeros serve as the quadrature nodes for the corresponding Gaussian
    /// rule. They are computed as the eigenvalues of the Jacobi tridiagonal
    /// matrix associated with this polynomial family.
    fn zeros(&self) -> Vec<F>;
}
