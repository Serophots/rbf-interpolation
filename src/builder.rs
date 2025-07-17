use std::iter::{repeat, Product};

use na::{
    ArrayStorage, Const, DefaultAllocator, DimMin, DimName, Matrix, RealField, Scalar,
    SquareMatrix, Storage, StorageMut, U1, Vector, allocator::Allocator,
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
    T,
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
    Multiquadratic { epsilon: T },
    /// 1/sqrt(1 + r^2)
    InverseMultiquadratic { epsilon: T },
    /// 1/(1 + r^2)
    InverseQuadratic { epsilon: T },
    /// exp(-r^2)
    Gaussian { epsilon: T },
}

impl<T, const DEGREE: usize, const MONOMIALS: usize, const POINTS: usize, const DIM: usize>
    RBFInterpolatorBuilder<T, DEGREE, MONOMIALS, POINTS, DIM>
where
    T: Scalar + RealField + Copy + Product,
{
    /// Constructs NxN matrix of phi(||a - b||) for each combination of points a, b
    fn construct_phi<SP>(
        &self,
        points: &Matrix<T, Const<DIM>, Const<POINTS>, SP>,
    ) -> SquareMatrix<
        T,
        Const<POINTS>,
        <DefaultAllocator as Allocator<Const<POINTS>, Const<POINTS>>>::Buffer<T>,
    >
    where
        SP: Storage<T, Const<DIM>, Const<POINTS>>,
    {
        let iter = points
            .column_iter()
            .flat_map(|point_a| points.column_iter().zip(repeat(point_a)))
            .map(|(point_a, point_b)| {
                self.kernel((point_a - point_b).map(|e| e.powi(2)).sum().sqrt())
            });

        SquareMatrix::<
            T,
            Const<POINTS>,
            <DefaultAllocator as Allocator<Const<POINTS>, Const<POINTS>>>::Buffer<T>,
        >::from_iterator(iter)
    }

    /// Add polynomial constraints to phi matrix (POINTS + MONOMIALS)x(POINTS + MONOMIALS)
    fn embelish_phi<SP, SH>(
        &self,
        points: &Matrix<T, Const<DIM>, Const<POINTS>, SP>,
        phi: SquareMatrix<T, Const<POINTS>, SH>,
    ) -> SquareMatrix<
        T,
        Const<{ POINTS + MONOMIALS }>,
        <DefaultAllocator as Allocator<
            Const<{ POINTS + MONOMIALS }>,
            Const<{ POINTS + MONOMIALS }>,
        >>::Buffer<T>,
    >
    where
        SP: Storage<T, Const<DIM>, Const<POINTS>>,
        SH: StorageMut<T, Const<POINTS>, Const<POINTS>>,
        Const<{ POINTS + MONOMIALS }>: DimName,
        DefaultAllocator: Allocator<Const<POINTS>, Const<POINTS>>
            + Allocator<Const<{ POINTS + MONOMIALS }>, Const<{ POINTS + MONOMIALS }>>,
    {
        let mut phi = phi.resize_generic(
            Const::<{ POINTS + MONOMIALS }>,
            Const::<{ POINTS + MONOMIALS }>,
            na::zero(),
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
                let value: T = exponent
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
        values: Vector<T, Const<POINTS>, SV>,
    ) -> Vector<
        T,
        Const<{ POINTS + MONOMIALS }>,
        <DefaultAllocator as Allocator<Const<{ POINTS + MONOMIALS }>, U1>>::Buffer<T>,
    >
    where
        SV: StorageMut<T, Const<POINTS>, U1>,
        Const<{ POINTS + MONOMIALS }>: DimName,
    {
        values.resize_generic(Const::<{ POINTS + MONOMIALS }>, U1, na::zero())
    }

    /// Solves (phi) * (weights) = (values), storing the solved weights into values
    /// Introduces new generic M to simplify the traits
    fn solve_phi<const M: usize, SH, SV>(
        &self,
        phi: SquareMatrix<T, Const<M>, SH>,
        values: &mut Vector<T, Const<M>, SV>,
    ) -> bool
    where
        SH: Storage<T, Const<M>, Const<M>>,
        SV: StorageMut<T, Const<M>, U1>,
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
        points: Matrix<T, Const<DIM>, Const<POINTS>, SP>,
        values: Vector<T, Const<POINTS>, SV>,
    ) -> Option<
        RBFInterpolator<
            T,
            DEGREE,
            MONOMIALS,
            POINTS,
            DIM,
            SP,
            ArrayStorage<T, { POINTS + MONOMIALS }, 1>,
        >,
    >
    where
        SP: Storage<T, Const<DIM>, Const<POINTS>>,
        SV: StorageMut<T, Const<POINTS>, U1>,
        //Convince compiler that phi is infact square
        Const<{ POINTS + MONOMIALS }>:
            DimMin<Const<{ POINTS + MONOMIALS }>, Output = Const<{ POINTS + MONOMIALS }>>,
    {
        let phi = self.construct_phi::<SP>(&points);
        let phi = self.embelish_phi::<SP, ArrayStorage<T, POINTS, POINTS>>(&points, phi);
        let mut values = self.embelish_values(values);
        if self.solve_phi(phi, &mut values) {
            Some(
                RBFInterpolator::<T, _, _, _, _, SP, ArrayStorage<T, { POINTS + MONOMIALS }, 1>> {
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
