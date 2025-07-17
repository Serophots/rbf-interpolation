use na::{Const, DimName, Matrix, Storage, U1, Vector};

use crate::{builder::RBFInterpolatorBuilder, powers::monomial_exponents};

type F = f64;

pub struct RBFInterpolator<
    const DEGREE: usize,
    const MONOMIALS: usize,
    const POINTS: usize,
    const DIM: usize,
    SP,
    SW,
> where
    SP: Storage<F, Const<DIM>, Const<POINTS>>,
    SW: Storage<F, Const<{ POINTS + MONOMIALS }>, U1>,
    Const<{ POINTS + MONOMIALS }>: DimName,
{
    pub(crate) kernel: RBFInterpolatorBuilder<DEGREE, MONOMIALS, POINTS, DIM>,
    pub(crate) points: Matrix<F, Const<DIM>, Const<POINTS>, SP>,
    pub(crate) weights: Vector<F, Const<{ POINTS + MONOMIALS }>, SW>,
}

impl<const DEGREE: usize, const MONOMIALS: usize, const POINTS: usize, const DIM: usize, SP, SW>
    RBFInterpolator<DEGREE, MONOMIALS, POINTS, DIM, SP, SW>
where
    SP: Storage<F, Const<DIM>, Const<POINTS>>,
    SW: Storage<F, Const<{ POINTS + MONOMIALS }>, U1>,
    Const<{ POINTS + MONOMIALS }>: DimName,
{
    pub fn interpolate<S1>(&self, point_b: &Vector<F, Const<DIM>, S1>) -> F
    where
        S1: Storage<F, Const<DIM>, U1>,
    {
        let mut weights = self.weights.iter();

        // Add the phi terms
        let phi: F = self
            .points
            .column_iter()
            .map(|point_a| {
                let weight = weights.next().unwrap();
                let phi = self
                    .kernel
                    .kernel((point_a - point_b).map(|e| e.powi(2)).sum().sqrt());

                weight * phi
            })
            .sum();

        // Add the polynomial terms
        let exponents = monomial_exponents::<DIM, DEGREE>();
        debug_assert_eq!(exponents.len(), MONOMIALS);

        let polynomial: F = exponents
            .iter()
            .map(|exponent| {
                debug_assert_eq!(exponent.len(), point_b.len());
                let weight = weights.next().unwrap();
                let value: F = exponent
                    .iter()
                    .zip(point_b.row_iter())
                    .map(|(&exponent, ordinate)| ordinate[(0, 0)].powi(exponent))
                    .product();

                weight * value
            })
            .sum();

        debug_assert!(weights.next().is_none());

        phi + polynomial
    }
}
