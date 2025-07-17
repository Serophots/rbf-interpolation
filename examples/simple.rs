use nalgebra::{Matrix, SMatrix, SVector, Vector2};
use rbf_interpolation::builder::RBFInterpolatorBuilder;

fn main() {
    //Define points (x,y)
    let points: SMatrix<f64, 2, 5> = Matrix::from_columns(&[
        Vector2::new(2.0, 2.0),
        Vector2::new(3.0, 4.0),
        Vector2::new(6.0, 4.0),
        Vector2::new(1.0, 1.0),
        Vector2::new(7.0, 7.0),
    ]);

    //Define values for each points (z)
    let values = SVector::<f64, 5>::new(2.0, 6.0, 4.0, 4.0, 5.0);

    //Construct the interpolant
    // Must provide the number of monomials terms for your choice of polynomial degree and points dimension
    // i.e.
    let interpolant = RBFInterpolatorBuilder::<f64, 1, 3, 5, 2>::ThinPlateSpline
        .build(points, values)
        .unwrap();

    interpolant.interpolate(&Vector2::new(10.0, 10.0));
}
