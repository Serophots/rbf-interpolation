use std::iter::{Product, Sum};

use na::{Const, DimName, Matrix, RealField, Scalar, Storage, U1, Vector};

use crate::{builder::RBFInterpolatorBuilder, powers::monomial_exponents};

pub struct RBFInterpolator<
    T,
    const DEGREE: usize,
    const MONOMIALS: usize,
    const POINTS: usize,
    const DIM: usize,
    SP,
    SW,
> where
    T: Scalar,
    SP: Storage<T, Const<DIM>, Const<POINTS>>,
    SW: Storage<T, Const<{ POINTS + MONOMIALS }>, U1>,
    Const<{ POINTS + MONOMIALS }>: DimName,
{
    pub(crate) kernel: RBFInterpolatorBuilder<T, DEGREE, MONOMIALS, POINTS, DIM>,
    pub(crate) points: Matrix<T, Const<DIM>, Const<POINTS>, SP>,
    pub(crate) weights: Vector<T, Const<{ POINTS + MONOMIALS }>, SW>,
}

impl<T, const DEGREE: usize, const MONOMIALS: usize, const POINTS: usize, const DIM: usize, SP, SW>
    RBFInterpolator<T, DEGREE, MONOMIALS, POINTS, DIM, SP, SW>
where
    T: Scalar + RealField + Copy + Sum + Copy + Product,
    SP: Storage<T, Const<DIM>, Const<POINTS>>,
    SW: Storage<T, Const<{ POINTS + MONOMIALS }>, U1>,
    Const<{ POINTS + MONOMIALS }>: DimName,
{
    pub fn interpolate<S1>(&self, point_b: &Vector<T, Const<DIM>, S1>) -> T
    where
        S1: Storage<T, Const<DIM>, U1>,
    {
        let mut weights = self.weights.iter();

        // Add the phi terms
        let phi: T = self
            .points
            .column_iter()
            .map(|point_a| {
                let &weight = weights.next().unwrap();
                let phi = self
                    .kernel
                    .kernel((point_a - point_b).map(|e| e.powi(2)).sum().sqrt());

                weight * phi
            })
            .sum();

        // Add the polynomial terms
        let exponents = monomial_exponents::<DIM, DEGREE>();
        debug_assert_eq!(exponents.len(), MONOMIALS);

        let polynomial: T = exponents
            .iter()
            .map(|exponent| {
                debug_assert_eq!(exponent.len(), point_b.len());
                let weight = weights.next().unwrap();
                let value: T = exponent
                    .iter()
                    .zip(point_b.row_iter())
                    .map(|(&exponent, ordinate)| ordinate[(0, 0)].powi(exponent))
                    .product();

                *weight * value
            })
            .sum();

        debug_assert!(weights.next().is_none());

        phi + polynomial
    }
}
