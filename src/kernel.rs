use crate::builder::RBFInterpolatorBuilder;

type F = f64;

impl<const DEGREE: usize, const MONOMIALS: usize, const N: usize, const D: usize>
    RBFInterpolatorBuilder<DEGREE, MONOMIALS, N, D>
{
    pub(crate) fn kernel(&self, r: F) -> F {
        match self {
            RBFInterpolatorBuilder::Linear => r,
            RBFInterpolatorBuilder::ThinPlateSpline => {
                if r == na::zero() {
                    na::zero()
                } else {
                    r.powi(2) * r.ln()
                }
            }
            RBFInterpolatorBuilder::Cubic => r.powi(3),
            RBFInterpolatorBuilder::Quintic => r.powi(5),
            RBFInterpolatorBuilder::Multiquadratic { .. } => todo!(),
            RBFInterpolatorBuilder::InverseMultiquadratic { epsilon, .. } => {
                na::one::<F>() / ((r / epsilon).powi(2) + na::one::<F>()).sqrt()
            }
            RBFInterpolatorBuilder::InverseQuadratic { .. } => todo!(),
            RBFInterpolatorBuilder::Gaussian { epsilon, .. } => {
                na::one::<F>() / ((r / epsilon).powi(2) + na::one::<F>()).exp()
            }
        }
    }
}
