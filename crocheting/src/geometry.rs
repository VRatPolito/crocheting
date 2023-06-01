pub mod curves {
    use nalgebra as na;
    use na::Unit;
    use nalgebra::{Point3, Vector3, UnitQuaternion};
    const NUM_SAMPLE_LUT: usize = 15;

    use crate::errors::errors::{CrochetError, ErrorType};

    use std::fs::File;
    use std::io::Write;

    use dyn_clone::DynClone;
    use dyn_clone;


    // CURVES
    // TRAIT ↓↓
    pub trait Curve: DynClone {
        fn evaluate(&self, u: f64) ->  Result<Point3<f64>, CrochetError>;
        fn first_derivative(&self, u: f64) -> Result<Vector3<f64>, CrochetError>;
        fn first_derivative_normalized(&self, u: f64) -> Result<Vector3<f64>, CrochetError>;
        fn second_derivative(&self, u: f64) -> Result<Vector3<f64>, CrochetError>;
        fn second_derivative_normalized(&self, u: f64) -> Result<Vector3<f64>, CrochetError>;

        //Versori che mi tornano comodi, non fanno torsioni
        fn versor_i(&self, u: f64) -> Result<Vector3<f64>, CrochetError>;
        fn versor_j(&self, u: f64) -> Result<Vector3<f64>, CrochetError>;
        fn versor_k(&self, u: f64) -> Result<Vector3<f64>, CrochetError>;


        // Versori dipendenti dalla direzione della derivata seconda. N B K dall'algebra
        fn t_tangent_versor(&self, u: f64) -> Result<Vector3<f64>, CrochetError>;
        fn n_normal_versor(&self, u: f64) -> Result<Vector3<f64>, CrochetError>;
        fn b_binormal_versor(&self, u: f64) -> Result<Vector3<f64>, CrochetError>;


        fn distance_from_start(&self, u: f64) -> f64; //Dato U, mi dice la distanza dall'inizio della curva
        fn distance_reversed(&self, dist: f64) ->  Result<f64, CrochetError>; //Quale U è a distanza "dist" dall'inizio?
        fn center_osculating_circle(&self, u: f64) -> Result<Point3<f64>, CrochetError>;//Quale è il centro della circ. osculatrice al punto u?
    
    
        // Metodo comodo per la visualizzazione
        fn to_file(&self, filename: Option<String>) -> String;
    }

    dyn_clone::clone_trait_object!(Curve);

    // Beziér
    #[derive(Clone)]
    pub struct Bezier {
        start: Point3<f64>,
        end: Point3<f64>,
        control_point_1: Option<Point3<f64>>,
        control_point_2: Option<Point3<f64>>,
        lut: Option<[f64; NUM_SAMPLE_LUT]>,
        lut_basis: Vec<(Vector3<f64>,Vector3<f64>)>
    }

    impl Bezier {
        pub fn new(
            start: Point3<f64>,
            end: Point3<f64>,
            control_point_1: Option<Point3<f64>>,
            control_point_2: Option<Point3<f64>>,
        ) -> Bezier {

            let mut bz = Bezier {
                start: start,
                end: end,
                control_point_1: control_point_1,
                control_point_2: control_point_2,
                lut: None,
                lut_basis: Vec::new()
            };
            
            bz.generate_lut();
            bz.generate_j_axis_lut();

            return bz;
        }

        /// Genero gli assi 
        fn generate_j_axis_lut(&mut self) -> ()
        {
            // Trovo un versore ortonormale tramite Gram-Schmidt
            let j_axis: Vector3<f64>;
            const X: Vector3<f64> = Vector3::new(1.0, 0.0, 0.0);
            const Y: Vector3<f64> = Vector3::new(1.0, 0.0, 0.0);
            let tang = self.versor_i(0.0).unwrap();
            if X.dot(&tang) > 0.001
            {
                j_axis = (X - X.dot(&tang)*tang).normalize();
            }
            else
            {
                j_axis = (Y - Y.dot(&tang)*tang).normalize();
            }

            
            let mut arr = vec![(self.first_derivative_normalized(0.0).unwrap(), j_axis)];

            for i in 1..NUM_SAMPLE_LUT
            {
                let curr_u: f64 = i as f64 / (NUM_SAMPLE_LUT as f64 -1.0);
                let past_u: f64 = (i-1) as f64 / (NUM_SAMPLE_LUT as f64 -1.0);
                let new_tan = self.versor_i(curr_u).unwrap();
                let old_tan = self.versor_i(past_u).unwrap();

                // let new_rot = UnitQuaternion::new(
                //     old_tan.cross(&new_tan)*std::f64::consts::PI / 2.0
                // );
                let new_rot: UnitQuaternion<f64> = UnitQuaternion::rotation_between_axis(
                    &Unit::new_unchecked(old_tan), &Unit::new_unchecked(new_tan)
                ).unwrap();

                let new_j: Vector3<f64> = new_rot * arr.last().unwrap().1;

                arr.push(
                    (
                        new_tan,
                        (new_j - new_j.dot(&new_tan)*new_tan).normalize()
                    )
                );
            }

            self.lut_basis = arr;
        }

        fn generate_lut(&mut self) -> () {
            if self.lut.is_none() {
                let mut new_lut = [0.0; NUM_SAMPLE_LUT];

                let mut distance = 0.0;
                let mut prev_point = Point3::new(self.start.x, self.start.y, self.start.z);
                for n in 0..NUM_SAMPLE_LUT {
                    let samplelut = NUM_SAMPLE_LUT as f64;
                    let v = self.evaluate(n as f64 / (samplelut - 1.0)).unwrap();
                    let new_point = Point3::new(v.x, v.y, v.z);
                    distance += na::distance(
                        &Point3::new(prev_point.x, prev_point.y, prev_point.z),
                        &new_point,
                    );
                    new_lut[n] = distance;

                    prev_point = new_point;
                }

                self.lut = Some(new_lut);
            }
        }
    }

    impl Curve for Bezier {
        fn evaluate(&self, u: f64) -> Result<Point3<f64>, CrochetError> {
            if u > 1.0
            {
                return Err(
                    CrochetError{
                        description: format!("Errore calcolo curva a u→ {}",u).to_string(),
                        error_type: ErrorType::OutOfRange
                    }
                )
            }
            let t = u.clamp(0.0, 1.0);

            if self.control_point_1.is_some()
                &&
                self.control_point_2.is_some()
            {
                let b1 = Bezier::new(
                    self.start,
                    self.control_point_2.unwrap(),
                    Some(self.control_point_1.unwrap()),
                    None
                );
                let b2 = Bezier::new(self.control_point_1.unwrap(),
                    self.end,
                    Some(self.control_point_2.unwrap()),
                    None
                );

                let res = add_points(
                    (1.0-t)*b1.evaluate(t)?,
                    t*b2.evaluate(t)?
                );

                return Ok(res);
            }
            else if self.control_point_1.is_some()
            {
                let res = add_points(
                    (1.0 - t) * add_points((1.0 - t) * self.start, t * self.control_point_1.unwrap()),
                    t * add_points((1.0 - t) * self.control_point_1.unwrap(), t * self.end)
                );                    
                return Ok(res);
            }
            else
            {
                let a = add_points( self.start * (1.0 - t), self.end * t);
                return Ok(Point3::new(a.x, a.y, a.z));
            }
        }


        fn versor_i(&self, u: f64) -> Result<Vector3<f64>, CrochetError> {
            return self.first_derivative_normalized(u);
        }
        fn versor_j(&self, u: f64) -> Result<Vector3<f64>, CrochetError> {
            if u > 1.0
            {
                return Err(
                    CrochetError{
                        description: "".to_string(),
                        error_type: ErrorType::OutOfRange
                    }
                )
            }

            let lut = self.lut_basis.clone();
            
            let scaled_u = u*(NUM_SAMPLE_LUT-1) as f64;
            let up_index = (scaled_u.ceil() as usize).clamp(0, NUM_SAMPLE_LUT-1);
            let low_index = (scaled_u.floor() as usize).clamp(0, NUM_SAMPLE_LUT-1);

            let factor = scaled_u - low_index as f64;

            let rot = UnitQuaternion::rotation_between_axis(
                &Unit::new_unchecked(
                    lut[low_index].0
                ),
                &Unit::new_unchecked(
                    lut[up_index].0
                )
            ).unwrap();

            return Ok(
                UnitQuaternion::identity().nlerp(&rot, factor) * lut[low_index].1
            );

            // FISSATO A ASSE X TEMPORANEAMENTE
            //return Ok(Vector3::new(1.0, 0.0, 0.0));
        }
        fn versor_k(&self, u: f64) -> Result<Vector3<f64>, CrochetError> {
            return Ok(
                self.versor_i(u)?.cross(&self.versor_j(u)?).normalize()
            );
        }


        // Algebra, T N e B di algebra, rivolti in base alla der. seconda.
        fn t_tangent_versor(&self, u: f64) -> Result<Vector3<f64>, CrochetError>
        {
            return self.first_derivative_normalized(u);
        }
        fn n_normal_versor(&self, mut u: f64) -> Result<Vector3<f64>, CrochetError>
        {
            if u > 1.0
            {
                u = 1.0;
            }
            let c = match self.center_osculating_circle(u)
            {
                Ok(el) => el,
                Err(el) => panic!("Panic per errore {}", el)
            };
            return Ok(
                    (c - self.evaluate(u).unwrap()).normalize()
                );
        }
        fn b_binormal_versor(&self, u: f64) -> Result<Vector3<f64>, CrochetError>
        {
            return Ok(
                self.t_tangent_versor(u)?.cross(&self.n_normal_versor(u)?).normalize()
            );
        }


        fn first_derivative(&self, mut u: f64) -> Result<Vector3<f64>, CrochetError>
        {
            const STEP: f64 = 0.05;
            if u > 1.0 - STEP {
                u -= STEP;
            }

            let p_0 = self.evaluate(u)?;
            let p_1 = self.evaluate(u + STEP)?;
            let res = (p_1 - p_0)/STEP;

            return Ok(res);
        }

        fn first_derivative_normalized(&self, u: f64) -> Result<Vector3<f64>, CrochetError> {
            Ok(self.first_derivative(u)?.normalize())
        }

        fn second_derivative(&self, mut u: f64) -> Result<Vector3<f64>, CrochetError>
        {
            const STEP: f64 = 0.05;
            if u > 1.0 - STEP {
                u -= STEP;
            }

            let p_0 = self.first_derivative(u)?;
            let p_1 = self.first_derivative(u + STEP)?;
            let res = (p_1 - p_0)/STEP;
            return Ok(res);
        }

        fn second_derivative_normalized(&self, u: f64) -> Result<Vector3<f64>, CrochetError> {
            let res = self.second_derivative(u)?;
            #[cfg(debug_assertions)]
            {
                println!("Versore î a t={}: {} \n", u, res);
            }
            

            Ok(res)
        }

        fn distance_from_start(&self, u: f64) -> f64
        {
            // Devo usare la LUT
            let my_lut = self.lut.unwrap();
            let index = u * (NUM_SAMPLE_LUT-1) as f64;
            
            let lower_index = index.floor();
            let upper_index = index.ceil();
            let fac = index - lower_index;

            return my_lut[lower_index as usize] * (1.0-fac) + my_lut[upper_index as usize]*fac;
        }

        fn distance_reversed(&self, dist: f64) -> Result<f64, CrochetError>
        {
            let lut = self.lut.unwrap();
            for i in 1..lut.len()
            {
                let elem = lut[i];
                let prev = lut[i-1];
                if dist >= prev && dist < elem
                {
                    let indice_calcolato = (i-1) as f64 + (dist - prev)/(elem - prev);
                    return Ok(indice_calcolato / (NUM_SAMPLE_LUT-1) as f64);
                }
            }
            // Dovrebbe restituire un errore o un risultato se è roba fuori scala
            return Err(CrochetError{description: "distance reversed out of range".to_string(), error_type:ErrorType::OutOfRange});
        }

        fn center_osculating_circle(&self, u: f64) -> Result<Point3<f64>, CrochetError>
        {
            // |r' X r''| / |r'|^3
            let first_der = self.first_derivative(u)?;
            let second_der = self.second_derivative(u)?;
            let k_curvature = (
                first_der.cross(&second_der).norm()
            )
            /
            (
                first_der.norm().powi(3)
            );

            if k_curvature < crate::crochet_consts::crochet_consts::ZERO_CURVATURE_TOLERANCE
            {
                return Err(CrochetError{description: "".to_string(), error_type: ErrorType::RadiusInfiniteDistance});
            }

            let radius = 1.0 / k_curvature;

            // P0 + second_der_versor*radius
            let center = self.evaluate(u)? + radius*second_der.normalize();

            return Ok(center);
        }

        fn to_file(&self, filename: Option<String>) -> String
        {
            let mut outp = "".to_string();
            for i in 0..=50
            {
                let p: Point3<f64> = self.evaluate(i as f64 / 50.0).unwrap();
                outp.push_str(&format!("{} {} {}\n", &p.x, &p.y, &p.z));
            }
            if filename.is_some()
            {
                let mut output = File::create(filename.unwrap()).unwrap();
                write!(output, "{}", outp).unwrap();
            }
            return outp;
        }
    }

    fn add_points(p1: Point3<f64>, p2: Point3<f64>) -> Point3<f64>
    {
        return Point3::new(
            p1.x + p2.x,
            p1.y + p2.y,
            p1.z + p2.z
        );
    }
}

pub mod shapes {
    use nalgebra::{Vector3, Rotation3, Translation3, Point3};
    use std::f64::consts::PI;
    use dyn_clone::DynClone;
    use crate::crochet_consts::crochet_consts::{get_width, get_v_rotation_surface};

    pub trait Shape: DynClone {
        fn ev(&self, v: f64) -> Point3<f64>;
        fn set_translation(&mut self, p: Point3<f64>) -> ();
        fn set_rotation_axis(&mut self, new_x: Vector3<f64>, new_y: Vector3<f64>, new_z: Vector3<f64>) -> ();
        fn v_max_shape_radius(&self) -> f64;

        /// Interpola le _shape_ con fattore tra 0 e 1, e alla coordinata _v_
        fn lerp(&self, factor: f64, v: f64, sh_2: &Box<dyn Shape>) -> Point3<f64>
        {
            let p_1 = self.ev(v);
            let p_2 = sh_2.ev(v);

            return Point3::new(
                p_1.x * (1.0-factor) + p_2.x * factor,
                p_1.y * (1.0-factor) + p_2.y * factor,
                p_1.z * (1.0-factor) + p_2.z * factor
            );
        }


        /// Distanza minima tra due Shape nello spazio
        fn min_distance(&self, sh_2: Box<dyn Shape>) -> f64{
            let mut min: f64 = f64::MAX;
            for i in 0..50
            {
                let p1= self.ev(i as f64 / 50.0);
                let p2= sh_2.ev(i as f64 / 50.0);
                
                let dist = nalgebra::distance(&p1, &p2);
                if dist < min
                {
                    min = dist;
                }
            }
            return min;
        }
    }

    dyn_clone::clone_trait_object!(Shape);


    // FLAT
    #[derive(Copy, Clone)]
    pub struct Circle {
        radius: f64,
        rotation: Rotation3<f64>,
        translation: Translation3<f64>
    }

    impl Circle {
        pub fn new(radius: f64) -> Circle {
            Circle {
                radius: radius,
                rotation: Rotation3::identity(),
                translation: Translation3::new(0.0, 0.0, 0.0)
            }
        }

        pub fn get_radius(&self) -> f64
        {
            self.radius
        }

        pub fn new_from_origin_and_versors(
            radius: f64,
            origin: Point3<f64>,
            i_j_k: &[Vector3<f64>; 3]
        ) -> Circle
        {
            let mut c = Circle::new(radius);
            c.set_translation(origin);
            c.set_rotation_axis(i_j_k[1], i_j_k[2], i_j_k[0]);

            return c;
        }
    }

    
    impl Shape for Circle {
        /// V è tra 0 e 1
        /// Restituisce le coordinate X e Y piane
        /// di dove passa la curva
        fn ev(&self, v: f64) -> Point3<f64> {
            if v < 0.0 || v > 1.0
            {
                panic!("Valore v: {}, fuori dall'intervallo 0-1", v);
            }
            let mut p = Point3::new((v * 2.0 * PI).cos()*self.radius, (v * 2.0 * PI).sin()*self.radius, 0.0);

            p = self.rotation * p;
            //println!("{}", p);
            p = Point3::new(
                self.translation.x + p.x,
                self.translation.y + p.y,
                self.translation.z + p.z
            );
            p
        }


        fn set_rotation_axis(&mut self, new_x: Vector3<f64>, new_y: Vector3<f64>, new_z: Vector3<f64>)
        {
            let x = new_x.normalize();
            let y = new_y.normalize();
            let z = new_z.normalize();
            self.rotation = Rotation3::from_basis_unchecked(
                &[
                x,
                y,
                z
            ]);
        }

        fn set_translation(&mut self, p: Point3<f64>)
        {
            self.translation = Translation3::new(p.x, p.y, p.z);
        }

        /// Quale è il raggio massimo?
        fn v_max_shape_radius(&self) -> f64
        {
            return self.radius;
            
            // Utile per shape più generale
            // let find_radius = |v: f64| {};
            // const ITER: u32 = 8;
            // let mut v_limits: (f64, f64) = (0.0, 1.0);
            // for i in 0..ITER
            // {
            //     let middle = (v_limits.0 + v_limits.1)/2.0;
                
            // }

            // todo!();
        }
    }

    

    #[derive(Clone)]
    pub struct MultiPointShape
    {
        points: Vec<(f64, Point3<f64>)>,
        rotation: Rotation3<f64>,
        translation: Translation3<f64>,
        lut_s: Vec<(f64, f64)>
    }
    impl MultiPointShape
    {
        /// Non inserire una seconda volta il valore a 1.0
        pub fn new(mut p: Vec<(f64, Point3<f64>)>) -> Self
        {
            p.push((1.0, p[0].1));

            let mut weights: Vec<f64> = Vec::new();
            for i in 0..p.len()-1
            {
                weights.push(nalgebra::distance(&p[i].1, &p[i+1].1));
            }
            let totale = weights.iter().copied().reduce(| a: f64, b: f64 | -> f64{ a+b }).unwrap();

            let mut new_lut: Vec<(f64, f64)> = Vec::new();

            new_lut.push((0.0, 0.0));

            {
                let mut accumulated: f64 = 0.0;
                for i in 0..weights.len()
                {
                    accumulated += weights[i];
                    new_lut.push(
                        (accumulated / totale, p[i+1].0)
                    );
                    
                }
            }
            

            MultiPointShape
            {
                points: p,
                rotation: Rotation3::identity(),
                translation: Translation3::new(0.0, 0.0, 0.0),
                lut_s: new_lut
            }
        }

        /// Costruisce una forma
        /// da un array di tuple
        /// 
        /// In ogni tupla (f64, f64), il primo valore è v tra 0 e 1,
        /// il secondo valore è il raggio in quel punto
        pub fn radial_shape(
            radius: Vec<(f64, f64)>
        ) -> Box<dyn Shape>
        {
            Box::new(MultiPointShape::new(
                radius.iter()
                .map(|e|{
                    (
                        e.0,
                        Point3::new(
                            (e.0 * 2.0 * PI).cos()*e.1, (e.0 * 2.0 * PI).sin()*e.1, 0.0
                        )
                    )
                }).collect()
            ))
        }

        pub fn new_ellipse(radius_a: f64, radius_b: f64, local_move: (f64, f64)) -> Self
        {
            const NUM_SUBDIVISIONS: usize = 40;
            let mut points: Vec<(f64, Point3<f64>)> = Vec::new();
            for i in 0..NUM_SUBDIVISIONS
            {
                let t = i as f64/NUM_SUBDIVISIONS as f64;
                let x = radius_a* (2.0*std::f64::consts::PI*t).cos() - local_move.0;
                let y = radius_b* (2.0*std::f64::consts::PI*t).sin() - local_move.1;

                points.push(
                    (t, Point3::new(x, y, 0.0))
                );
            }
            return MultiPointShape::new(points);
        }

        pub fn new_magic_ring(translation: (f64, f64) ) -> Self
        {
            let w = get_width();

            let circonferenza = 6.0*w;

            let radius = circonferenza / (2.0*PI);

            let c = MultiPointShape::new_ellipse(radius, radius, translation);
            return c;
        }

        fn lut_compute(&self, v: f64) -> f64
        {
            for i in 0..self.lut_s.len()-1
            {
                if v >= self.lut_s[i].0 && v <= self.lut_s[i+1].0
                {
                    let factor = (v - self.lut_s[i].0) / (self.lut_s[i+1].0 - self.lut_s[i].0);
                    return self.lut_s[i].1*(1.0-factor) + self.lut_s[i+1].1*factor;
                }
            }
            panic!("Errore lut multipointshape");
        }
    }
    impl Shape for MultiPointShape
    {
        fn ev(&self, v_old: f64) -> Point3<f64>
        {
            let v_diff = get_v_rotation_surface();

            let v: f64;
            {
                let mut tmp_v = v_old - v_diff;
                if tmp_v < 0.0
                {
                    tmp_v += 1.0;
                }
                
                tmp_v = tmp_v - tmp_v.floor();

                v=tmp_v;
            }

            let v = self.lut_compute(v);
            // Dovrebbe usare parametrizzazione intrinseca
            let mut p = self.points[self.points.len()-1].1;
            //let index_limits: (usize, usize) = (0, 1);
            for i in 0..self.lut_s.len()-1
            {
                if v >= self.points[i].0 && v <= self.points[i+1].0
                {
                    let factor = (v - self.points[i].0) / (self.points[i+1].0 - self.points[i].0);
                    p = Point3::new(
                        self.points[i].1.x*(1.0 - factor) + self.points[i+1].1.x*(factor),
                        self.points[i].1.y*(1.0 - factor) + self.points[i+1].1.y*(factor),
                        self.points[i].1.z*(1.0 - factor) + self.points[i+1].1.z*(factor)
                    );
                }
            }
            p = self.rotation * p;
            p = Point3::new(
                self.translation.x + p.x,
                self.translation.y + p.y,
                self.translation.z + p.z
            );
            return p;
        }



        fn set_translation(&mut self, p: Point3<f64>)
        {
            self.translation = Translation3::new(p.x, p.y, p.z);
        }

        fn set_rotation_axis(&mut self, new_x: Vector3<f64>, new_y: Vector3<f64>, new_z: Vector3<f64>)
        {
            let x = new_x.normalize();
            let y = new_y.normalize();
            let z = new_z.normalize();
            self.rotation = Rotation3::from_basis_unchecked(
                &[
                x,
                y,
                z
            ]);
        }


        /// Quale è il raggio massimo?
        /// 
        /// Restituisce raggio
        fn v_max_shape_radius(&self) -> f64
        {
            // Utile per shape più generale
            let mut max = Vector3::new(
                self.translation.x - self.points[0].1.x,
                self.translation.y - self.points[0].1.y,
                self.translation.z - self.points[0].1.z
            ).norm();
            
            for i in 1..self.points.len()
            {
                let n = Vector3::new(
                    self.translation.x - self.points[i].1.x,
                    self.translation.y - self.points[i].1.y,
                    self.translation.z - self.points[i].1.z
                ).norm();
                if n > max
                {
                    max = n;
                }
            }
            
            return max;
        }
    }
}

pub mod transform
{
    use nalgebra::{Point3};
    use crate::param_surface::surface::{UVSurface, TauCoords};
    use crate::crochet_consts;

    #[derive(Debug)]
    pub struct TriLinearInterpolation
    {
        point_0_0_0: Point3<f64>,
        point_x_0_0: Point3<f64>,
        point_0_y_0: Point3<f64>,
        point_x_y_0: Point3<f64>,
        point_0_0_z: Point3<f64>,
        point_x_0_z: Point3<f64>,
        point_0_y_z: Point3<f64>,
        point_x_y_z: Point3<f64>,
    }
    impl TriLinearInterpolation
    {
        pub fn new_on_surface(points: [TauCoords; 4], surf: &Box<dyn UVSurface>) -> Self
        {
            let spessore = crochet_consts::crochet_consts::get_thickness();
            let below_face: Vec<Point3<f64>> = points.iter().map(
                |c| {
                    surf.to_3d_point(*c).unwrap()
                    -
                    surf.normal_versor(*c)*spessore
                }
            ).collect();
            let upper_face: Vec<Point3<f64>> = points.iter().map(
                |c| {
                    surf.to_3d_point(*c).unwrap()
                    +
                    surf.normal_versor(*c)*spessore
                }
            ).collect();

            return TriLinearInterpolation::new(below_face.try_into().unwrap(), upper_face.try_into().unwrap());
        }

        /// Restituisce la faccia sotto e faccia sopra
        pub fn get_faces(&self) -> (Vec<Point3<f64>>,Vec<Point3<f64>>)
        {
            let below_face: Vec<Point3<f64>> = vec![
                self.point_0_0_0,
                self.point_x_0_0,
                self.point_x_y_0,
                self.point_0_y_0
            ];
            let up_face: Vec<Point3<f64>> = vec![
                self.point_0_0_z,
                self.point_x_0_z,
                self.point_x_y_z,
                self.point_0_y_z
            ];
            return (
                below_face,
                up_face
            );
        }

        /// Restituisce il punto di controllo derivata corretto rispetto
        /// alla faccia indicata
        pub fn correct_control_point(&self, punto_controllo_derivata: Point3<f64>, punto: Point3<f64>) -> Point3<f64>
        {
            let trasformato: Point3<f64> = self.transform_point(punto);
            let derivata_trasformata: Point3<f64> = self.transform_point(punto_controllo_derivata);

            let versore_normale = match punto {
                punto if punto.x < 0.001 => {
                    (self.point_0_y_0 - self.point_0_0_0).cross(
                        &(self.point_0_0_z - self.point_0_0_0)
                    ).normalize()
                },
                punto if punto.y < 0.001 => {
                    - (self.point_x_0_0 - self.point_0_0_0).cross(
                        &(self.point_0_0_z - self.point_0_0_0)
                    ).normalize()
                },
                punto if punto.x > 0.99 => {
                    - (self.point_0_y_0 - self.point_0_0_0).cross(
                        &(self.point_0_0_z - self.point_0_0_0)
                    ).normalize()
                },
                punto if punto.y > 0.99 => {
                    (self.point_x_0_0 - self.point_0_0_0).cross(
                        &(self.point_0_0_z - self.point_0_0_0)
                    ).normalize()
                },
                _ => {
                    return derivata_trasformata;
                }
            };
            let dir_v = derivata_trasformata - trasformato;

            return 
                trasformato
                +
                dir_v.dot(&versore_normale)*versore_normale;
        }

        /// Passare le coordinate delle facce in senso antiorario
        pub fn new(below_face: [Point3<f64>;4], up_face: [Point3<f64>;4]) -> Self
        {
            return TriLinearInterpolation{
                point_0_0_0: below_face[0],
                point_x_0_0: below_face[1],
                point_0_y_0: below_face[3],
                point_x_y_0: below_face[2],
                point_0_0_z: up_face[0],
                point_x_0_z: up_face[1],
                point_0_y_z: up_face[3],
                point_x_y_z: up_face[2],
            }
        }

        pub fn transform(&self, u: f64, v: f64, w: f64) -> Point3<f64>
        {
            fn lerp(p1: Point3<f64>, p2: Point3<f64>, factor: f64) -> Point3<f64>
            {
                return Point3::new(
                    p1.x*(1.0-factor) + p2.x*factor,
                    p1.y*(1.0-factor) + p2.y*factor,
                    p1.z*(1.0-factor) + p2.z*factor
                );
            }
            let plane_down: Point3<f64> = lerp(
                lerp(self.point_0_0_0, self.point_0_y_0, v),
                lerp(self.point_x_0_0, self.point_x_y_0, v),
                u
            );
            let plane_up: Point3<f64> = lerp(
                lerp(self.point_0_0_z, self.point_0_y_z, v),
                lerp(self.point_x_0_z, self.point_x_y_z, v),
                u
            );

            return lerp(
                plane_down,
                plane_up,
                w
            );
        }

        pub fn transform_point(&self, p: Point3<f64>) -> Point3<f64>
        {
            return self.transform(p.x, p.y, p.z);
        }
    }
}