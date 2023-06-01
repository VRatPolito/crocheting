#[macro_use]
extern crate lazy_static;

pub mod nodes;
pub mod param_surface;
pub mod geometry;
pub mod crochet_consts;
mod errors;
pub mod instructions;
pub mod mesh_viz;


#[cfg(test)]
mod tests {
    use crate::geometry::transform::TriLinearInterpolation;
    use nalgebra::{Point3, Vector3};
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }

    #[test]
    fn trilinearInterp()
    {
        let tr = TriLinearInterpolation::new(
            [
                Point3::new(1.0, 1.0, 1.0),
                Point3::new(1.0, 3.0, 1.0),
                Point3::new(0.0, 3.0, 2.0),
                Point3::new(0.0, 1.0, 2.0)
            ],
            [
                Point3::new(1.0, 1.0, 2.0),
                Point3::new(1.0, 3.0, 2.0),
                Point3::new(0.0, 3.0, 3.0),
                Point3::new(0.0, 1.0, 3.0)
            ]
        );
        assert!(
            check_close(
                tr.transform(0.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 1.0)
            )
        );
        assert!(
            check_close(
                tr.transform(1.0, 1.0, 1.0),
                Point3::new(0.0, 3.0, 3.0)
            )
        );
        assert!(
            check_close(
                tr.transform(0.5, 0.5, 0.0),
                Point3::new(0.5, 2.0, 1.5)
            )
        );
        assert!(
            check_close(
                tr.transform(0.5, 0.5, 1.0),
                Point3::new(0.5, 2.0, 2.5)
            )
        );

        assert!(
            check_close(
                tr.transform(0.0, 0.0, 0.5),
                Point3::new(1.0, 1.0, 1.5)
            )
        );
        
        //tr.transform(, v: f64, w: f64)
    }

    #[test]
    fn trilinearInterpDerivata()
    {
        let tr = TriLinearInterpolation::new(
            [
                Point3::new(4.0, -1.0, 0.0),
                Point3::new(4.0, 1.0, 0.0),
                Point3::new(4.0, 4.0, 3.0),
                Point3::new(4.0, 0.0, 3.0),
            ],
            [
                Point3::new(2.0, -1.0, 0.0),
                Point3::new(2.0, 1.0, 0.0),
                Point3::new(0.0, 4.0, 3.0),
                Point3::new(0.0, 0.0, 3.0),
            ]
        );


        let p = tr.correct_control_point(
            Point3::new(0.5, 0.5, 0.5 )
            , Point3::new(0.5, 0.0, 0.5));
        let vec_to_check = p - tr.transform_point(Point3::new(0.5, 0.0, 0.5));
        assert_eq!(
            vec_to_check.dot(&Vector3::new(1.0, 0.0, 0.0)),
            0.0
        );
        assert_eq!(
            vec_to_check.dot(&Vector3::new(0.0, 1.0, 0.0)),
            0.0
        );
    }

    fn check_close(p1: Point3<f64>, p2: Point3<f64>) -> bool
    {
        return nalgebra::distance(&p1, &p2) < 0.0000001;
    }

}