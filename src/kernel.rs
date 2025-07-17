use na::{RealField, Scalar};

use crate::builder::RBFInterpolatorBuilder;

impl<T, const DEGREE: usize, const MONOMIALS: usize, const N: usize, const D: usize>
    RBFInterpolatorBuilder<T, DEGREE, MONOMIALS, N, D>
where
    T: Scalar + RealField + Copy,
{
    pub(crate) fn kernel(&self, r: T) -> T {
        match self {
            RBFInterpolatorBuilder::Linear => -r,
            RBFInterpolatorBuilder::ThinPlateSpline => {
                if r == na::zero() {
                    na::zero()
                } else {
                    r.powi(2) * r.ln()
                }
            }
            RBFInterpolatorBuilder::Cubic => r.powi(3),
            RBFInterpolatorBuilder::Quintic => -r.powi(5),
            RBFInterpolatorBuilder::Multiquadratic { epsilon } => {
                -((r * *epsilon).powi(2) + na::one::<T>()).sqrt()
            }
            RBFInterpolatorBuilder::InverseMultiquadratic { epsilon, .. } => {
                na::one::<T>() / ((r * *epsilon).powi(2) + na::one::<T>()).sqrt()
            }
            RBFInterpolatorBuilder::InverseQuadratic { epsilon } => {
                na::one::<T>() / ((r * *epsilon).powi(2) + na::one::<T>())
            }
            RBFInterpolatorBuilder::Gaussian { epsilon, .. } => {
                na::one::<T>() / ((r * *epsilon).powi(2) + na::one::<T>()).exp()
            }
        }
    }
}
