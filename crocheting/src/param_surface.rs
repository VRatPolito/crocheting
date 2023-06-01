//  Crocheting - a Rust library to computer knitting instructions on parametric surfaces
//  Copyright © 2022 Massimo Gismondi
//
//  This program is free software: you can redistribute it and/or
//  modify it under the terms of the GNU General Public License
//  as published by the Free Software Foundation, either version 3
//  of the License, or (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//  GNU General Public License for more details.

//  You should have received a copy of the GNU General Public License
//  along with this program. If not, see <https://www.gnu.org/licenses/>. 

pub mod surface
{
    use nalgebra::{Point3, Vector3};
    use nalgebra as na;
    use crate::geometry::curves::{Curve};
    use crate::geometry::shapes::{Circle, Shape, MultiPointShape};
    use crate::crochet_consts::crochet_consts;
    use crate::errors::errors::{CrochetError, ErrorType};
    use std::fs::File;
    use std::io::Write;
    use std::ops::{Add, Sub};
    use dyn_clone::DynClone;
    use dyn_clone;
    
    // TRAITS
    pub trait UVSurface: DynClone
    {
        /// Scorrendo tau, trovo tutti i punti dell'elica
        fn helix(&self, tau: f64) -> TauCoords
        {
            TauCoords {
                tau: tau,
                v: tau - tau.floor()
            }
        }

        fn to_triangles(&self, starting_index: i32) -> (Vec<i32>, Vec<Point3<f64>>)
        {
            // let subdivision_tau: usize = 4;
            let subdivision_tau: usize = (
                (
                    self.get_curve().distance_from_start(1.0) / crochet_consts::get_yarn_radius()
                ).ceil() / 4.0
            ) as usize + 1;
            const SUBDIVISION_V: usize = 6;

            let mut indices: Vec<i32> = Vec::new();
            let mut points: Vec<Point3<f64>> = Vec::new();

            // Genero prima i vertici
            for i in 0..=subdivision_tau
            {
                let tau = (i as f64 / subdivision_tau as f64)*self.tau_max_value();
                for j in 0..SUBDIVISION_V
                {
                    let v = j as f64 / SUBDIVISION_V as f64;

                    points.push(
                        self.to_3d_point(TauCoords{tau: tau, v: v}).unwrap()
                    )
                }
            }

            fn compute_index(tau: i32, v: i32, max_v: usize) -> i32
            {
                max_v as i32 * tau + v
            }

            for i in 0..subdivision_tau
            {
                let tau_low: i32 = i as i32;
                let tau_up: i32 = (i+1) as i32;
                for j in 0..SUBDIVISION_V
                {
                    let v_prev: i32 = j as i32;
                    let v_next: i32 = match j == SUBDIVISION_V-1
                    {
                        true => {0},
                        false => {j as i32 + 1}
                    } ;

                    
                    let quad = [
                        compute_index(tau_low, v_prev, SUBDIVISION_V),
                        compute_index(tau_low, v_next, SUBDIVISION_V),
                        compute_index(tau_up, v_next, SUBDIVISION_V),
                        compute_index(tau_up, v_prev, SUBDIVISION_V)
                    ];

                    // Primo triangolo  
                    indices.extend(
                        [
                            quad[0],
                            quad[3],
                            quad[1]
                        ]
                    );

                    // Secondo triangolo  
                    indices.extend(
                        [
                            quad[1],
                            quad[3],
                            quad[2]
                        ]
                    );
                }
            }

            for i in 0..indices.len()
            {
                indices[i] = indices[i] + starting_index;
            }

            return (indices, points);
        }

        /// Converte coordinate UV nelle X/Y/Z corrispondenti
        /// Utile a scopo di visualizzazione su grafici
        /// oppure per misurare distanze con pitagora
        fn to_3d_point(&self, coord: TauCoords) -> Result<Point3<f64>, CrochetError>;
    
        /// Nuovo punto a una maglia di distanza in altezza
        /// o "dist"
        fn height_increment(&self, coord: TauCoords, dist: f64) -> Result<TauCoords, CrochetError>
        {
            let u = self.tau_to_u(coord.tau)?;
            // 
            // La superficie può allargarsi o stringersi nel salire.
            // Devo tenerne conto
            // Tramite ds/dtau
            let d_tau: f64 = 0.001;
            let p_base = coord;
            let p_aum_tau = coord + TauCoords{tau: d_tau, v: 0.0};
            let p_aum_v = coord + TauCoords{tau: 0.0, v: d_tau};
            let p_aum_helix = coord + TauCoords{tau: d_tau, v: d_tau};
            let geodesic_dist_tau = self.geodesic_dist(p_base, p_aum_tau, None)?;
            let geodesic_dist_v = self.geodesic_dist(p_base, p_aum_v, None)?;

            let direzione = self.get_curve().versor_i(u)?.normalize();
            let versore_tau_global: Vector3<f64> = (self.to_3d_point(p_aum_tau)? - self.to_3d_point(p_base)?).normalize();
            let versore_v_global: Vector3<f64> = (self.to_3d_point(p_aum_v)? - self.to_3d_point(p_base)?).normalize();
            let versore_v_spirale: Vector3<f64> = (self.to_3d_point(p_aum_helix)? - self.to_3d_point(p_base)?).normalize();

            let normal: Vector3<f64> = self.normal_versor(coord);

            let orto_to_spiral_global: Vector3<f64> = - versore_v_spirale.cross(&normal);

            
            
            // DI quanto tau devo spostarmi per ottenere
            // la lunghezza dist voluta?
            // ds_su_dtau è una derivata che andrò a mettere a denominatore
            // Ogni tau, quanto mi cambia la distanza lungo l'elica?
            let ds_su_dtau = geodesic_dist_tau / d_tau;
            let ds_su_v = geodesic_dist_v / d_tau;

            let new_coord: TauCoords;
            if crochet_consts::get_diagonal_height_increment()
            {
                // NUOVO
                fn square_interp(v: f64, max_v: (f64, f64)) -> (f64, f64)
                {
                    assert!(v>=0.0 && v<=1.0);
                    let limits = (max_v.1, -max_v.1);
                    if v<0.25
                    {
                        return (v*4.0*max_v.0, max_v.1);
                    }
                    else if v<0.75
                    {
                        // println!("limiti v {:?}", limits);
                        let v = (v-0.25)*2.0;
                        return (max_v.0, (1.0 - v)*limits.0 + v*limits.1);
                    }
                    else
                    {
                        let v = (v-0.75)*4.0;
                        return ((1.0-v)*max_v.0, -max_v.1);
                    }
                }
                let max_delta_coordinates = (dist / ds_su_dtau, -dist / ds_su_v);
                const ITERATIONS: usize = 10;
                let mut bin_search_limits: (f64, f64) = (0.0, 1.0);

                loop
                {
                    let middle_factor = (bin_search_limits.0 + bin_search_limits.1)/2.0;
                    let interp = square_interp(middle_factor, max_delta_coordinates);
                    let middle_coords = coord + TauCoords{tau: interp.0, v: interp.1};
                    let middle_direction: Vector3<f64> = (self.to_3d_point(middle_coords)? - self.to_3d_point(coord)?).normalize();

                    let proj = middle_direction.dot(&versore_v_spirale);
                    if proj < 0.0
                    {
                        bin_search_limits = (middle_factor, bin_search_limits.1);
                    }
                    else
                    {
                        bin_search_limits = (bin_search_limits.0, middle_factor);
                    }

                    if proj.abs() < 0.01
                    {
                        break;
                    }
                }

                
                let found_factor = (bin_search_limits.0 + bin_search_limits.1)/2.0;


                let found_delta_coords = square_interp(found_factor, max_delta_coordinates);
                let target_coords_not_scaled = coord+TauCoords{tau: found_delta_coords.0, v: found_delta_coords.1};
                let scaled_delta_coords = coord.lerp(
                    &target_coords_not_scaled,
                    (crochet_consts::get_height() / self.geodesic_dist(coord, target_coords_not_scaled, None)?).clamp(0.0, 1.0)
                );
                new_coord = scaled_delta_coords;
            // FINE NUOVO
            }
            else
            {
                // VECCHIO
                let direct_tau_incremented_coord = TauCoords{
                    tau: coord.tau + dist / ds_su_dtau,
                    v: coord.v
                };
                let spostamento = (self.to_3d_point(direct_tau_incremented_coord)? - self.to_3d_point(coord)?);

                let projection_on_v = (spostamento).dot(&versore_v_global);

                new_coord = direct_tau_incremented_coord + TauCoords{tau: 0.0, v: - projection_on_v / ds_su_v};
                // FINE vecchio
            }
            
            

            
            
            // Fa ritornare a un errore se siamo fuori dalla superficie
            self.tau_to_u(new_coord.tau)?;
            //println!("{:?} {:?}", coord, new_coord);
            return Ok(new_coord);

            return Ok(
                TauCoords{
                    tau: coord.tau + dist / ds_su_dtau,
                    v: coord.v
                }
            );
        }

        fn get_curve(&self) -> &Box<dyn Curve>;

        fn get_tau_lut(&self) -> &Vec<f64>;
        
        fn normal_versor(&self, t: TauCoords) -> Vector3<f64>
        {
            const DELTA: f64 = 0.001;
            let start_coord: TauCoords;
            if t.tau+DELTA>self.tau_max_value()
            {
                start_coord = TauCoords{tau: t.tau - DELTA, v:t.v} ;
            }
            else
            {
                start_coord = t;
            }
            
            let origine = self.to_3d_point(start_coord).unwrap();

            return (self.to_3d_point(start_coord+TauCoords{tau:0.0, v: DELTA}).unwrap() - origine).normalize().cross(
                &(self.to_3d_point(start_coord+TauCoords{tau:DELTA, v: 0.0}).unwrap() - origine).normalize()
            ).normalize();
        }

        /// Tau è un remap di U
        fn tau_to_u(&self, tau: f64) -> Result<f64, CrochetError>
        {
            //let tau_scaled = tau * (self.tau_lut.len()as f64-1.0) ;
            if tau > self.tau_max_value() as f64
            {
                return Err(
                    CrochetError{
                        description: format!("value {} out of bounds", tau).to_string(),
                        error_type: ErrorType::OutOfRange
                    }
                );
            }
            let u_low = self.get_tau_lut()[tau.floor() as usize];
            let u_up = self.get_tau_lut()[tau.ceil() as usize];
            let p = self.helix(tau).v;
            let u = (u_low*(1.0-p) + u_up*p).clamp(0.0, 1.0);

            return Ok(u);
        }

        /// Remap inverso, utile in alcuni conti quando bisogna tornare
        /// da U indietro verso tau
        fn u_to_tau(&self, u: f64) -> f64
        {
            let lut = self.get_tau_lut();
            for i in 1..lut.len()
            {
                let elem = lut[i];
                let prev = lut[i-1];
                if u >= prev && u <= elem
                {
                    let indice_calcolato = (i-1) as f64 + (u - prev)/(elem - prev);
                    return indice_calcolato;
                }
            }
            return lut[lut.len()-1];
        }

        /// Può servire sapere da fuori quanto è il massimo
        /// valore di tau, per non eccederlo
        fn tau_max_value(&self) -> f64
        {
            return self.u_to_tau(1.0);
        }

        /// Utile a scopo di visualizzazione
        fn to_file(&self, filename: Option<String>) -> String
        {
            let mut outp = "".to_string();
            let points = self.tau_max_value().floor() as i32*20;
            for i in 1..points
            {
                let step: f64 = self.tau_max_value() / points as f64;
                let tau: f64 = (i as f64)*step;
                let p = match self.to_3d_point(self.helix(tau))
                {
                    Ok(e) => e,
                    _ => break
                };
                outp.push_str(&format!("{} {} {}\n", &p.x, &p.y, &p.z));
            }
            if filename.is_some()
            {
                let mut output = File::create(filename.unwrap()).unwrap();
                write!(output, "{}", outp).unwrap();
            }
            return outp;
        }



        /// Distanza geodesica, mi muovo lungo la curva
        /// Viene approssimato con sample_n campioni
        /// 
        /// Se sample_n non è dato, viene usato un valore
        /// predefinito abbastanza basso
        fn geodesic_dist(&self, from: TauCoords, to: TauCoords, sample_n: Option<u32>) -> Result<f64, CrochetError>
        {
            let sample = sample_n.unwrap_or(3);
            let mut dist: f64 = 0.0;
            for i in 1..=sample
            {
                let prev = from.lerp(&to, (i as f64-1.0) / sample as f64);
                let current =  from.lerp(&to, (i as f64) / sample as f64);
                dist += na::distance(
                    &self.to_3d_point(prev)?,
                    &self.to_3d_point(current)?
                );
            }

            Ok(dist)
        }

        /// Nuovo punto a una maglia di distanza in larghezza
        fn width_increment(&self, start_coord: TauCoords, dist: f64) -> Result<TauCoords, CrochetError>
        {
            // Tramite ds/dtau
            let d_tau: f64 = 0.01;
            let p1 = start_coord;
            let p2 = start_coord + TauCoords{tau: d_tau, v: d_tau};
            let geodesic_dist = self.geodesic_dist(p1, p2, None);
            
            // DI quanto tau devo spostarmi per ottenere
            // la lunghezza dist voluta?
            // ds_su_dtau è una derivata che andrò a mettere a denominatore
            // Ogni tau, quanto mi cambia la distanza lungo l'elica?
            let ds_su_dtau = geodesic_dist? / d_tau;


            return Ok(
                start_coord + TauCoords{
                tau: dist / ds_su_dtau,
                v: dist / ds_su_dtau
            }) ;
        }
    }
    dyn_clone::clone_trait_object!(UVSurface);
    
    #[derive(Debug)]
    #[derive(Copy, Clone, PartialEq)]
    pub struct TauCoords
    {
        pub tau: f64,// valore intero a ogni inizio ciclo. Rimappatura di U
        pub v: f64// valore di angolo tra 0 e 1
    }
    impl TauCoords
    {
        pub fn lerp(&self, other: &TauCoords, factor: f64) -> TauCoords
        {
            if factor < 0.0 || factor > 1.0
            {
                panic!("Factor {} out of bounds 0-1", factor);
            }

            let new_tau = self.tau*(1.0-factor) + other.tau*factor;
            let mut v1 = self.v;
            let mut v2 = other.v;
            if (v2-v1).abs() > 0.5
            {
                if v2 > v1
                {
                    v1 += 1.0;
                }
                else
                {
                    v2 += 1.0;
                }
            }

            let new_v = v1*(1.0-factor) + v2*factor;

            TauCoords{
                tau: new_tau,
                v: new_v - new_v.floor()
            }
        }

        /// Delta_v corrisponde a una semplice sottrazione
        /// tra le "v" di due taucoords.
        /// 
        /// La funzione si occupa di gestire risultati
        /// minori di -0.5 o maggiori di 0.5 
        /// che derivano dalla sottrazione diretta.
        /// 
        /// Attenzione, usare solo per piccole variazioni di delta_v
        /// 
        /// ```
        /// use crocheting::param_surface::surface::TauCoords;
        /// let a = TauCoords{tau: 0.8, v: 0.9};
        /// let b = TauCoords{tau: 0.85, v: 0.1};
        /// 
        /// let coord_to_move = TauCoords{tau: 0.95, v: 0.95};
        /// 
        /// let delta_v = b.v - a.v;
        /// let res: TauCoords = coord_to_move.move_right_helix(delta_v);
        /// 
        /// println!("{:?}, {}", res, delta_v);
        /// 
        /// assert_eq!(
        ///     res.tau,
        ///     1.15
        /// );
        /// assert!(
        ///     (res.v - 0.15).abs() < 0.00000001
        /// );
        ///```
        pub fn move_right_helix(&self, delta_v: f64) -> TauCoords
        {
            let patched_delta_v: f64;
            if delta_v < -0.5
            {
                patched_delta_v = 1.0 + delta_v;
            }
            else if delta_v > 0.5
            {
                patched_delta_v = delta_v - 1.0;
            }
            else
            {
                patched_delta_v = delta_v;
            }

            return self.clone() + TauCoords{tau: patched_delta_v, v: patched_delta_v};
        }
    }

    impl Add for TauCoords
    {
        type Output = Self;
        fn add(self, other: Self) -> Self {
            let mut new_v = self.v + other.v;
            if new_v < 0.0
            {
                new_v += 1.0;
            }
            Self {
                tau: self.tau + other.tau,
                v: new_v - new_v.floor(),
            }
        }
    }

    impl Sub for TauCoords
    {
        type Output = Self;
        fn sub(self, other: Self) -> Self {
            let mut new_v = self.v - other.v;
            if new_v < 0.0
            {
                new_v += 1.0;
            }
            Self {
                tau: self.tau - other.tau,
                v: new_v - new_v.floor(),
            }
        }
    }

    // Implemento i traits
    #[derive(Clone)]
    pub struct Cylinder
    {
        r1: f64,
        r2: f64,
        curve: Box<dyn Curve>,
        tau_lut: Vec<f64>,//Associazione variabile tau → u
    }
    
    impl Cylinder
    {
        fn min_circle_dist(c1: &Circle, c2: &Circle) -> f64
        {
            let mut dist = -1.0;
            const SAMPLE: i32 = 20;
            for i in 0..SAMPLE
            {
                let p1 = c1.ev(
                    i as f64 / (SAMPLE + 1) as f64
                );
                let p2 = c2.ev(
                    i as f64 / (SAMPLE + 1) as f64
                );
                let d = na::distance(&p1, &p2);
                
                if dist < 0.0 || d < dist
                {
                    dist = d;
                }
            }
            
            return dist;
        }

        pub fn new(curve: Box<dyn Curve>, radius1: f64, radius2: f64) -> Cylinder
        {
            let mut new_lut: Vec<f64> = Vec::new();
            new_lut.push(0.0);


            let mut t: f64 = 0.0;
            while t < 1.0
            {
                // Devo inizializzare la LUT di tau→u

                let mut prev_shape = Circle::new(radius1*(1.0-t) + radius2*t);
                prev_shape.set_translation(curve.evaluate(t).unwrap());
                prev_shape.set_rotation_axis(curve.versor_j(t).unwrap(), curve.versor_k(t).unwrap(), curve.versor_i(t).unwrap());
                // Creo il cerchio successivo
                let mut circles_dist = 0.0;
                while circles_dist < crochet_consts::get_height()
                {
                    if t>1.0
                    {
                        break;
                    }
                    // 'sta roba è una porcata
                    let cur_radius = radius1*(1.0-t) + radius2*t;
                    let mut cur_circle = Circle::new(cur_radius);
                    cur_circle.set_translation(curve.evaluate(t).unwrap());
                    cur_circle.set_rotation_axis(curve.versor_j(t).unwrap(), curve.versor_k(t).unwrap(), curve.versor_i(t).unwrap());

                    circles_dist = Cylinder::min_circle_dist(&prev_shape, &cur_circle);
                    
                    t+=0.0001;
                }
                
                #[cfg(debug_assertions)]
                {
                    // for i in 0..20
                    // {
                    //     //let p1 = prev_shape.ev(i as f64 / (20 + 1) as f64);
                    //     //println!("{} {} {}", p1.x, p1.y, p1.z);
                    // }
                }
                

                new_lut.push(t);
            }

            // L'ultimo valore sarà sopra a 1 sicuro... lo rimuovo
            //new_lut.pop();

            //println!("{:?}", new_lut);


            Cylinder{
                r1: radius1,
                r2: radius2,
                curve: curve,
                tau_lut: new_lut,
                //helix_distance_lut: None
            }
        }

        pub fn new_no_taulut(curve: Box<dyn Curve>, radius1: f64, radius2: f64) -> Cylinder
        {
            return Cylinder
            {
                curve: curve,
                r1: radius1,
                r2: radius2,
                tau_lut: vec![0.0, 1.0]
            };
        }
    }

    impl UVSurface for Cylinder
    {
        fn get_curve<'a>(&'a self) -> &'a Box<dyn Curve>
        {
            &self.curve
        }

        fn get_tau_lut<'a>(&'a self) -> &'a Vec<f64>
        {
            &self.tau_lut
        }

        fn to_3d_point(&self, coord: TauCoords) -> Result<Point3<f64>, CrochetError>
        {
            let u = self.tau_to_u(coord.tau)?;

            // Middle circle
            let origin = self.curve.evaluate(u)?;
            let ijk = [
                self.curve.versor_i(u)?,
                self.curve.versor_j(u)?,
                self.curve.versor_k(u)?
                ];
            
            let cur_circle = Circle::new_from_origin_and_versors(
                self.r1*(1.0-u) + self.r2*u,
                origin,
                &ijk
            );
            
            return Ok(cur_circle.ev(coord.v));
        }

        

    }

    /// Superficie che riceve un numero arbitrario di _shape_ 
    /// e di _u_ corrispondenti sulla curva
    #[derive(Clone)]
    pub struct MultiShapeSurface
    {
        curve: Box<dyn Curve>,
        shapes: Vec<(f64, Box<dyn Shape>)>,
        tau_lut: Vec<f64>
    }
    impl MultiShapeSurface
    {
        pub fn new(
            curve: Box<dyn Curve>,
            mut start_shapes: Vec<(f64, Box<dyn Shape>)>
        ) -> Self
        {
            

            // Sposto ogni shape lungo la curva, ruotandola opportunamente
            for e in start_shapes.iter_mut()
            {
                e.1.set_translation( curve.evaluate( e.0 ).unwrap() );
                e.1.set_rotation_axis(
                    curve.versor_j(e.0).unwrap(),
                    curve.versor_k(e.0).unwrap(),
                    curve.versor_i(e.0).unwrap()
                );
            }
            
            let mut t: f64 = 0.0;
            let mut new_lut: Vec<f64> = vec![t];

            while t < 1.0
            {
                // Trovo un punto sopra al valore richiesto
                let prev: Box<dyn Shape> = TmpSurface::shape_at(&start_shapes, t);
                let mut incremented_t = t;
                while incremented_t < 1.0
                {
                    if prev.min_distance( TmpSurface::shape_at(&start_shapes, incremented_t)) > crochet_consts::get_height()
                    {
                        break;
                    }
                    incremented_t += 0.005;
                }
                if incremented_t >= 1.0
                {
                    break;
                }

                // Binary search per avvicinarmici velocemente
                let mut limits: (f64, f64) = (t, incremented_t);
                const ITERATIONS: i32 = 8;
                for _i in 0..ITERATIONS
                {
                    let middle = (limits.0 + limits.1)/2.0;
                    if prev.min_distance( TmpSurface::shape_at(&start_shapes, middle)) > crochet_consts::get_height()
                    {
                        limits.1 = middle;
                    }
                    else
                    {
                        limits.0 = middle;
                    }
                }

                new_lut.push(limits.1);
                t = limits.1;
            }


            // Continuo linearmente
            let last = new_lut[new_lut.len()-1];
            let prev_last = new_lut[new_lut.len()-2];
            let dist = last-prev_last;

            let new_t = last+dist;

            new_lut.push(match new_t >= 1.0
            {
                true => new_t,
                false => 1.0
            });


            MultiShapeSurface{
                curve: curve,
                shapes: start_shapes,
                tau_lut: new_lut
            }
        }
    }


    
    impl UVSurface for MultiShapeSurface
    {
        fn get_curve(&self) -> &Box<dyn Curve>
        {
            &self.curve
        }

        fn get_tau_lut(&self) -> &Vec<f64>
        {
            &self.tau_lut
        }

        fn to_3d_point(&self, coord: TauCoords) -> Result<Point3<f64>, CrochetError>
        {
            let u = self.tau_to_u(coord.tau)?;         
            return Ok(TmpSurface::evaluate(&self.shapes, u, coord.v));
        }
    }



    /// Superficie indicizzata con u e v anziché tau
    /// Utile per inizializzare tau_lut
    struct TmpSurface{}
    impl TmpSurface
    {
        pub fn evaluate(shapes: &Vec<(f64, Box<dyn Shape>)>, u: f64, v: f64) -> Point3<f64>
        {
            let s = TmpSurface::find_shapes_and_factor(shapes, u);

            return s.0.lerp(s.2, v, s.1);
        }

        pub fn shape_at(shapes: &Vec<(f64, Box<dyn Shape>)>, u: f64) -> Box<dyn Shape>
        {
            let mut new_shape_points: Vec<(f64, Point3<f64>)> = Vec::new();

            for i in 0..=40
            {
                let v = i as f64 / 40.0;
                new_shape_points.push(
                    (v, TmpSurface::evaluate(shapes, u, v))
                );
            }

            return Box::new(MultiPointShape::new(new_shape_points));
        }

        fn find_shapes_and_factor(shapes: &Vec<(f64, Box<dyn Shape>)>, u: f64) -> (&Box<dyn Shape>, &Box<dyn Shape>, f64)
        {
            for i in 0..shapes.len()-1
            {
                if u >= shapes[i].0 && u <= shapes[i+1].0
                {
                    let factor = (u - shapes[i].0) / (shapes[i+1].0 - shapes[i].0);
                    return (
                        &shapes[i].1,
                        &shapes[i+1].1,
                        factor
                    );
                }
            }
            
            panic!("Valore {} fuori scala", u);
        }
    }
}
