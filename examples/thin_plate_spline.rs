use gnuplot::{AutoOption::Auto, AxesCommon, Figure, HELIX};
use nalgebra::{Matrix, SMatrix, SVector, Vector2};
use rbf_interpolation::builder::RBFInterpolatorBuilder;

fn main() {
    let points: SMatrix<f64, 2, 5> = Matrix::from_columns(&[
        Vector2::new(2.0, 2.0),
        Vector2::new(3.0, 4.0),
        Vector2::new(6.0, 4.0),
        Vector2::new(1.0, 1.0),
        Vector2::new(7.0, 7.0),
    ]);
    let values = SVector::<f64, 5>::new(2.0, 6.0, 4.0, 4.0, 5.0);

    let interpolant = RBFInterpolatorBuilder::<1, 3, 5, 2>::ThinPlateSpline
        .build(points.clone(), values.clone())
        .unwrap();

    //15x15 with 0.1 resolution
    let mut surface = Vec::with_capacity(150*150);
    for i in 0..150 {
        for j in 0..150 {
            surface.push(interpolant.interpolate(&Vector2::new(i as f64 / 10.0, j as f64 / 10.0)));
        }
    }

    let mut points = points.row_iter();

    let mut fg = Figure::new();
    fg.axes3d()
        .set_title("Thin plate spline, 1 degree polynomial", &[])
        .surface(surface, 150, 150, Some((0.0, 0.0, 15.0, 15.0)), &[])
        .points(points.next().unwrap(), points.next().unwrap(), values.iter(), &[])
        .set_x_label("X", &[])
        .set_y_label("Y", &[])
        .set_z_label("Z", &[])
        .set_z_range(Auto, Auto)
        .set_palette(HELIX)
        .set_view(45.0, 175.0);
    fg.show().unwrap();
}
