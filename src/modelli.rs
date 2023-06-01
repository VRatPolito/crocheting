pub mod modelli
{
    use crocheting::param_surface::surface::{UVSurface, Cylinder, MultiShapeSurface};
    use crocheting::geometry::curves::Bezier;
    use crocheting::geometry::shapes::{MultiPointShape,Shape, Circle};
    use nalgebra::Point3;
    pub fn beastie() -> Box<dyn UVSurface>
    {
        let b = Bezier::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 7.0),
            None,
            None
        );
        
        let forme: Vec<(f64, Box<dyn Shape>)> = vec![
            (0.0,   Box::new(MultiPointShape::new_ellipse(3.0, 3.0, (0.0, 0.0))) ),
            (0.2,   Box::new(MultiPointShape::new_ellipse(3.0, 3.0, (0.0, 0.0))) ),
            (0.5,   Box::new(MultiPointShape::new_ellipse(3.0, 3.0, (3.0, 0.0))) ),
            (0.8,   Box::new(MultiPointShape::new_ellipse(3.0, 3.0, (0.0, 0.0))) ),
            (1.0,   Box::new(MultiPointShape::new_ellipse(3.0, 3.0, (0.0, 0.0))) )
        ];

        let s = MultiShapeSurface::new(Box::new(b), forme);

        return Box::new(s);
    }
}