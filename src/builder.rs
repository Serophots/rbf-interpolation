use std::{iter::repeat, usize};

use na::{
    ArrayStorage, Const, DefaultAllocator, DimMin, DimName, Matrix, SquareMatrix, Storage,
    StorageMut, U1, Vector, allocator::Allocator,
};

use crate::{powers::monomial_exponents, rbf::RBFInterpolator};

/// You must provide the number of monomials terms
/// for your choice of polynomial degree and points
/// dimension.
///
/// i.e. A degree 1 polynomial in two variables will
/// incur the following monomials: 1, x, y
/// Degree 2 will incur: 1, x, y, xy, x^2, y^2
///
/// In general, monomials = (degree + dimension) choose (degree)
pub enum RBFInterpolatorBuilder<
    const DEGREE: usize,
    const MONOMIALS: usize,
    const POINTS: usize,
    const DIM: usize,
> {
    /// r
    Linear,
    /// r^2 * log(r)
    ThinPlateSpline,
    /// r^3
    Cubic,
    /// r^5
    Quintic,
    /// -sqrt(1 + r^2)
    Multiquadratic { epsilon: F },
    /// 1/sqrt(1 + r^2)
    InverseMultiquadratic { epsilon: F },
    /// 1/(1 + r^2)
    InverseQuadratic { epsilon: F },
    /// exp(-r^2)
    Gaussian { epsilon: F },
}

type F = f64;

impl<const DEGREE: usize, const MONOMIALS: usize, const POINTS: usize, const DIM: usize>
    RBFInterpolatorBuilder<DEGREE, MONOMIALS, POINTS, DIM>
{
    /// Constructs NxN matrix of phi(||a - b||) for each combination of points a, b
    fn construct_phi<SP>(
        &self,
        points: &Matrix<F, Const<DIM>, Const<POINTS>, SP>,
    ) -> SquareMatrix<
        F,
        Const<POINTS>,
        <DefaultAllocator as Allocator<Const<POINTS>, Const<POINTS>>>::Buffer<F>,
    >
    where
        SP: Storage<F, Const<DIM>, Const<POINTS>>,
    {
        let iter = points
            .column_iter()
            .flat_map(|point_a| points.column_iter().zip(repeat(point_a)))
            .map(|(point_a, point_b)| {
                self.kernel((point_a - point_b).map(|e| e.powi(2)).sum().sqrt())
            });

        SquareMatrix::<
            F,
            Const<POINTS>,
            <DefaultAllocator as Allocator<Const<POINTS>, Const<POINTS>>>::Buffer<F>,
        >::from_iterator(iter)
    }

    /// Add polynomial constraints to phi matrix (POINTS + MONOMIALS)x(POINTS + MONOMIALS)
    fn embelish_phi<SP, SH>(
        &self,
        points: &Matrix<F, Const<DIM>, Const<POINTS>, SP>,
        phi: SquareMatrix<F, Const<POINTS>, SH>,
    ) -> SquareMatrix<
        F,
        Const<{ POINTS + MONOMIALS }>,
        <DefaultAllocator as Allocator<
            Const<{ POINTS + MONOMIALS }>,
            Const<{ POINTS + MONOMIALS }>,
        >>::Buffer<F>,
    >
    where
        SP: Storage<F, Const<DIM>, Const<POINTS>>,
        SH: StorageMut<F, Const<POINTS>, Const<POINTS>>,
        Const<{ POINTS + MONOMIALS }>: DimName,
        DefaultAllocator: Allocator<Const<POINTS>, Const<POINTS>>
            + Allocator<Const<{ POINTS + MONOMIALS }>, Const<{ POINTS + MONOMIALS }>>,
    {
        let mut phi = phi.resize_generic(
            Const::<{ POINTS + MONOMIALS }>,
            Const::<{ POINTS + MONOMIALS }>,
            0.0,
        );
        let exponents = monomial_exponents::<DIM, DEGREE>();

        assert_eq!(
            exponents.len(),
            MONOMIALS,
            "The choices of generics polynomial DEGREE and number of variables DIM did not match the number of MONONOMIALS."
        );

        for (i, point) in points.column_iter().enumerate() {
            for (j, exponent) in exponents.iter().enumerate() {
                debug_assert_eq!(exponent.len(), point.len());
                let value: F = exponent
                    .iter()
                    .zip(point.row_iter())
                    .map(|(&exponent, ordinate)| ordinate[(0, 0)].powi(exponent))
                    .product();

                phi[(i, POINTS + j)] = value;
                phi[(POINTS + j, i)] = value;
            }
        }

        phi
    }

    /// Add polynomial constraints to the values vector (POINTS + MONOMIALS)
    fn embelish_values<SV>(
        &self,
        values: Vector<F, Const<POINTS>, SV>,
    ) -> Vector<
        F,
        Const<{ POINTS + MONOMIALS }>,
        <DefaultAllocator as Allocator<Const<{ POINTS + MONOMIALS }>, U1>>::Buffer<F>,
    >
    where
        SV: StorageMut<F, Const<POINTS>, U1>,
        Const<{ POINTS + MONOMIALS }>: DimName,
    {
        values.resize_generic(Const::<{ POINTS + MONOMIALS }>, U1, 0.0)
    }

    /// Solves (phi) * (weights) = (values), storing the solved weights into values
    /// Introduces new generic M to simplify the traits
    fn solve_phi<const M: usize, SH, SV>(
        &self,
        phi: SquareMatrix<F, Const<M>, SH>,
        values: &mut Vector<F, Const<M>, SV>,
    ) -> bool
    where
        SH: Storage<F, Const<M>, Const<M>>,
        SV: StorageMut<F, Const<M>, U1>,
        Const<M>: DimMin<Const<M>, Output = Const<M>>,
        DefaultAllocator: Allocator<<Const<M> as DimMin<Const<M>>>::Output>,
    {
        phi.lu().solve_mut(values)
    }

    /// Number of points with dimension provided as a matrix of
    /// collum vectors, with their values in a seperate vector.
    ///
    /// Returns None when the linear system was not solveable.
    /// Will panic when the choice of added polynomial DEGREE
    /// and number of corresponding MONOMIAL terms are incompatible.
    /// Should satisfy: MONOMIAL = (DIM+DEGREE) choose DEGREE.
    pub fn build<SP, SV>(
        self,
        points: Matrix<F, Const<DIM>, Const<POINTS>, SP>,
        values: Vector<F, Const<POINTS>, SV>,
    ) -> Option<
        RBFInterpolator<
            DEGREE,
            MONOMIALS,
            POINTS,
            DIM,
            SP,
            ArrayStorage<F, { POINTS + MONOMIALS }, 1>,
        >,
    >
    where
        SP: Storage<F, Const<DIM>, Const<POINTS>>,
        SV: StorageMut<F, Const<POINTS>, U1>,
        //Convince compiler that phi is infact square
        Const<{ POINTS + MONOMIALS }>:
            DimMin<Const<{ POINTS + MONOMIALS }>, Output = Const<{ POINTS + MONOMIALS }>>,
    {
        let phi = self.construct_phi::<SP>(&points);
        let phi = self.embelish_phi::<SP, ArrayStorage<F, POINTS, POINTS>>(&points, phi);
        let mut values = self.embelish_values(values);
        if self.solve_phi(phi, &mut values) {
            Some(
                RBFInterpolator::<_, _, _, _, SP, ArrayStorage<F, { POINTS + MONOMIALS }, 1>> {
                    kernel: self,
                    points,
                    weights: values,
                },
            )
        } else {
            None
        }
    }
}
