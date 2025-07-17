#[cfg(test)]
mod tests {
    use std::time::Instant;

    use na::{Matrix, SMatrix, SVector, Vector2};

    use crate::builder::RBFInterpolatorBuilder;

    #[test]
    fn test() {
        let instant = Instant::now();

        let points: SMatrix<f64, 2, 5> = Matrix::from_columns(&[
            Vector2::new(2.0, 2.0),
            Vector2::new(3.0, -4.0),
            Vector2::new(6.0, -4.0),
            Vector2::new(-1.0, 1.0),
            Vector2::new(7.0, 7.0),
        ]);
        let values = SVector::<f64, 5>::new(2.0, 6.0, 4.0, 4.0, 5.0);

        let interpolant = RBFInterpolatorBuilder::<1, 3, 5, 2>::ThinPlateSpline
            .build(points, values)
            .unwrap();

        println!("construction {:?}", instant.elapsed());
        let instant = Instant::now();

        println!("test {}", interpolant.interpolate(&Vector2::new(2.0, 2.0)));
        println!("test {}", interpolant.interpolate(&Vector2::new(3.0, -4.0)));

        println!("interpolate {:?}", instant.elapsed());
    }
}
