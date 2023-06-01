use crocheting::param_surface::surface::{Cylinder, UVSurface, MultiShapeSurface};
use crocheting::geometry::curves::{Bezier, Curve};
use crocheting::geometry::shapes::{MultiPointShape, Circle, Shape};
use crocheting::mesh_viz::mesh_viz::StitchMesh;
use nalgebra::{Point3};
use crocheting::nodes::nodes;
use crocheting::nodes::nodes::{KnittableGraph, EdgeType};
use crocheting::instructions::instructions::InstructionsGenerator;

use std::time::{Instant, Duration};
use std::fs;

mod modelli;

fn main()
{
    std::fs::create_dir_all("dati_grafici").unwrap();
    unsafe
    {
        crocheting::crochet_consts::crochet_consts::set_stitch_dimensions(0.01, 0.01);
    }
    
    let b =Bezier::new(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.5, 0.0),
        None,
        None
    );
    #[cfg(debug_assertions)]
    {
        b.to_file(Some("dati_grafici/build_curve.mat".to_string()));
    }
    
    
    let s = Cylinder::new(Box::new(b), 0.05, 0.02);

    // let b =Bezier::new(
    //     Point3::new(0.0, 0.0, 0.0),
    //     Point3::new(0.0, 0.1, 0.0),
    //     Some(Point3::new(0.0, 0.05, 0.0)),
    //     //None,
    //     None
    // );
    // b.to_file(Some("curve.mat".to_string()));
    
    // let s = Cylinder::new(Box::new(b), 0.08, 0.03);
    #[cfg(debug_assertions)]
    {
        s.to_file(Some("dati_grafici/build_spirale_target.mat".to_string()));
    }

    let mut gr = nodes::CrochetGraph2::new(Box::new(s));
    gr.generate_graph();
    

    #[cfg(debug_assertions)]
    {
        gr.dot_graph(Some("dati_grafici/graph.dot".to_string()));

        gr.to_string_vertices(Some("dati_grafici/build_vertici.mat".to_string()));

        // let mut istr = InstructionsGenerator::new(Box::new(gr));
        // istr.human_instructions(Some("dati_grafici/build_istruzioni_testuali_generate.txt".to_string()));
    }
    
}


////////////////////////////////
// Sezione dei test
#[cfg(test)]
mod tests
{
    use std::fs::File;
    use std::io::Write;
    use crocheting::param_surface::surface::{Cylinder, TauCoords, UVSurface};
    #[test]
    fn knitting_test()
    {
        // println!("----------\n geodesic");
        // geodesic();
        // println!("----------\n cilindro_semplice");
        // cilindro_semplice();
        println!("----------\n knitted_test_1");
        knitted_test_1();
        // println!("----------\n knitted_test_2");
        // knitted_test_2();
        // println!("----------\n knitted_test_2_no_shortrow");
        // knitted_test_2_no_shortrow();
        // println!("----------\n bezier_examples");
        // bezier_examples();
        // println!("----------\n knitted_cylinder");
        // knitted_cylinder();
        // println!("----------\n godot_meshviz");
        // godot_meshviz();
        // println!("----------\n capunaman");
        // capunaman();
        // println!("----------\n capunaman_knitted_test_2");
        // capunaman_knitted_test_2();
        // println!("----------\n capunaman_multishape");
        // println!("capunaman_multishape crochetgraph");
        // multishape_1_crochetgraph();
        // println!("capunaman_multishape capunaman");
        // multishape_1_capunaman();
        // println!("----------\n oca celeste");
        // oca_celeste();
        // println!("----------\n beastie");
        // beastie();
        // println!("----------\n Spiral check");
        // spiral_check();

        // plots_for_compliance();


        // Time benchmark
        // println!("----------\n time_benchmark");
        // time_benchmark_godot_meshviz();
        // time_benchmark_oca_celeste();
        mean_variance_over_stitch_size();
    }

    fn dist_similar(n: f64, l: f64, precision: f64) -> bool
    {
        (n-l).abs() < precision
    }

    use super::*;

    #[test]
    fn test_straight_curve()
    {
        let b =Bezier::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 0.1),
            None,
            //None,
            None
        );
        assert_eq!(b.distance_from_start(1.0), 0.1);
    }

    #[test]
    fn test_curve_distance()
    {
        let b =Bezier::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 0.1),
            None,
            //None,
            None
        );
        for i in 0..20
        {
            let t: f64 = i as f64/20.0;
            assert!(
                dist_similar(t, b.distance_reversed(b.distance_from_start(t)).unwrap(), 0.0001)
            );
        }
    }

    #[test]
    fn test_tau_lut()
    {
        let b =Bezier::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 0.1),
            None,
            //None,
            None
        );
        let s = Cylinder::new(Box::new(b), 0.05, 0.10);
        let points = s.tau_max_value().floor() as i32*2;
        for i in 0..points
        {
            let step: f64 = s.tau_max_value() / points as f64;
            let t: f64 = i as f64*step;
            let u_calcolata = s.tau_to_u(t).unwrap();
            let tau_calcolata = s.u_to_tau(u_calcolata);
            assert!(
                dist_similar(
                    t,
                    tau_calcolata,
                    0.0001),
                    "Errore per tau={}, massimo={}\n{}",t,s.tau_max_value(),s.u_to_tau(s.tau_to_u(t).unwrap())
            );
        }
    }

    #[test]
    fn test_straight_surface()
    {
        let b =Bezier::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 0.1),
            None,
            //None,
            None
        );
        let s = Cylinder::new(Box::new(b), 0.05, 0.05);
    }

    
    fn test_cylinder_render()
    {
        std::fs::create_dir_all("dati_grafici/test_cylinder_render").unwrap();
        let b =Bezier::new(
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
            Some(Point3::new(0.0, 0.0, 0.0)),
            //None,
            Some(Point3::new(0.0, 0.0, 0.0)),
        );
        b.to_file(Some("dati_grafici/test_cylinder_render/spirale.mat".into()));
        let s = Cylinder::new(Box::new(b), 0.05, 0.05);

        let t = s.to_triangles(0);

        let mut output = File::create("dati_grafici/test_cylinder_render/out.mat").unwrap();
        for v in t.1.iter()
        {
            write!(output, "{} {} {}\n", v.x, v.y, v.z).unwrap();
        }
    }

    fn geodesic()
    {
        unsafe
        {
            crocheting::crochet_consts::crochet_consts::set_stitch_dimensions(0.05, 0.05);
        }
        let h = crocheting::crochet_consts::crochet_consts::get_height();
        let b =Bezier::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, h*5.0),
            None,
            //None,
            None
        );
        let s = Cylinder::new(Box::new(b), 0.05, 0.05);
        assert!(
            dist_similar(
                s.geodesic_dist(
                    TauCoords{tau: 1.0, v: 0.1},
                    TauCoords{tau: 1.0, v: 0.9}
                    , None).unwrap(), 
                s.geodesic_dist(
                    TauCoords{tau: 1.0, v: 0.9},
                    TauCoords{tau: 1.0, v: 0.1}
                    , None
                ).unwrap(),
                0.00001
            )
        );
    }


    /// Cilindro semplice
    fn cilindro_semplice()
    {
        unsafe
        {
            crocheting::crochet_consts::crochet_consts::set_stitch_dimensions(0.05, 0.05);
        }
        std::fs::create_dir_all("dati_grafici/cilindro_semplice").unwrap();
        let b =Bezier::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 0.5),
            // Some(Point3::new(0.0, 0.0, 0.45)),
            None,
            None
        );

        b.to_file(Some("dati_grafici/cilindro_semplice/curve.mat".to_string()));
        
        
        let s = Cylinder::new(Box::new(b), 0.1, 0.1);
    

        s.to_file(Some("dati_grafici/cilindro_semplice/spirale_target.mat".to_string()));

    
        let mut gr = nodes::CrochetGraph2::new(Box::new(s));
        gr.generate_graph();

        gr.to_string_vertices(Some("dati_grafici/cilindro_semplice/vertici.mat".to_string()));

        gr.start_rows_nodes(Some("dati_grafici/cilindro_semplice/startrows.mat".to_string()));
    }


    /// Questo è il primo test che cucito all'uncinetto come esempio
    /// in lana arancione
    fn knitted_test_1()
    {
        const SCALE: f64 = 16.0;
        unsafe
        {
            crocheting::crochet_consts::crochet_consts::set_stitch_dimensions(0.2, 0.2);
        }
        std::fs::create_dir_all("dati_grafici/knitted_test_1/").unwrap();
        let b =Bezier::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.5, 0.0)*SCALE,
            Some(Point3::new(0.0, 0.0, 0.45)*SCALE),
            //None,
            None
        );
        b.to_file(Some("dati_grafici/knitted_test_1/curve.mat".to_string()));
        
        let s = Cylinder::new(Box::new(b), 0.05*SCALE, 0.02*SCALE);
        s.to_file(Some("dati_grafici/knitted_test_1/spirale_target.mat".to_string()));

        {
            let mut gr: Box<dyn KnittableGraph> = Box::new(nodes::CrochetGraph2::new(Box::new(s.clone())));
            gr.generate_graph();
            gr.dot_graph(Some("dati_grafici/knitted_test_1/graph.dot".to_string()));
            gr.to_string_vertices(Some("dati_grafici/knitted_test_1/vertici.mat".to_string()));
            gr.start_rows_nodes(Some("dati_grafici/knitted_test_1/startrows.mat".to_string()));
            gr.edge_error_distribution(Some("dati_grafici/knitted_test_1/error_distribution_gismondi.mat".into()), None);
            gr.metrics_skew(Some("dati_grafici/knitted_test_1/skew_distribution_gismondi.mat".into()));
    
            let mut st = StitchMesh::new_from_graph(&gr);
            st.array_verts(false);
            st.data_for_blender_viz("dati_grafici/knitted_test_1/blender_facedata.txt");
    
            let istr = InstructionsGenerator::new(&st);
            istr.text_instructions(Some("dati_grafici/knitted_test_1/istruzioni.txt"));
        }
        {
            // Capunaman
            let mut gr_cap: Box<dyn KnittableGraph> = Box::new(nodes::CapunamanPaperGraph::new(Box::new(s.clone())));
            gr_cap.generate_graph();
            gr_cap.edge_error_distribution(Some("dati_grafici/knitted_test_1/error_distribution_capunaman.mat".into()), None);
            gr_cap.metrics_skew(Some("dati_grafici/knitted_test_1/skew_distribution_capunaman.mat".into()));
            let mut st = StitchMesh::new_from_graph(&gr_cap);
            st.array_verts(false);
            st.data_for_blender_viz("dati_grafici/knitted_test_1/blender_facedata_capunaman.txt");
    
            let istr = InstructionsGenerator::new(&st);
            istr.text_instructions(Some("dati_grafici/knitted_test_1/istruzioni_capunaman.txt"));
        }


        // static_run_benchmark(Box::new(s), "knitted_test_1", false);
    }


    /// Questo è il primo test che cucito all'uncinetto come esempio
    /// in lana arancione
    fn knitted_test_2_no_shortrow()
    {
        unsafe
        {
            crocheting::crochet_consts::crochet_consts::set_stitch_dimensions(0.02, 0.02);
            crocheting::crochet_consts::crochet_consts::set_flag_generate_short_row_flag(false);
            crocheting::crochet_consts::crochet_consts::set_snap_treshold(0.0);
        }
        std::fs::create_dir_all("dati_grafici/knitted_test_2_no_shortrow").unwrap();


        let b =Bezier::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.5, 0.0),
            Some(Point3::new(0.0, 0.0, 0.45)),
            //None,
            None
        );

        b.to_file(Some("dati_grafici/knitted_test_2_no_shortrow/curve.mat".to_string()));
        
        
        let s = Cylinder::new(Box::new(b), 0.05, 0.02);
    

        s.to_file(Some("dati_grafici/knitted_test_2_no_shortrow/spirale_target.mat".to_string()));

        // Coordinate anelli sovrapposti
        const SAMPLE_ANELLO: u32 = 30;
        let mut my_string = "".to_string();
        for i in 0..s.tau_max_value().floor() as i32
        {
            for n in 0..=30
            {
                let p = s.to_3d_point(
                    TauCoords{tau: i as f64, v: n as f64 / SAMPLE_ANELLO as f64}
                ).unwrap();
                my_string.push_str(&format!("{} {} {}\n", &p.x, &p.y, &p.z));
            }
        }
        
        let mut output = File::create("dati_grafici/knitted_test_2_no_shortrow/anelli_sovrapposti.mat").unwrap();
        write!(output, "{}", my_string).unwrap();

        

    
        let mut gr = nodes::CrochetGraph2::new(Box::new(s));
        gr.generate_graph();
        
        gr.start_rows_nodes(Some("dati_grafici/knitted_test_2_no_shortrow/startrows.mat".to_string()));

        gr.dot_graph(Some("dati_grafici/knitted_test_2_no_shortrow/graph.dot".to_string()));

        gr.to_string_vertices(Some("dati_grafici/knitted_test_2_no_shortrow/vertici.mat".to_string()));

        gr.debug_shortrow_level_nodes(Some("dati_grafici/knitted_test_2_no_shortrow/vertici_shortrow.mat".to_string()), 1);

        gr.edge_error_distribution(Some("dati_grafici/knitted_test_2_no_shortrow/error_distribution.mat".into()), None);
        gr.metrics_skew(Some("dati_grafici/knitted_test_2_no_shortrow/skew_distribution.mat".into()));

        unsafe
        {
            crocheting::crochet_consts::crochet_consts::set_flag_generate_short_row_flag(true);
            crocheting::crochet_consts::crochet_consts::reset_snap_treshold();
        }
    }


    /// Test 2, è il primo ma con una minore risoluzione.
    /// Esporta anche la versione ad anelli sovrapposti
    fn knitted_test_2()
    {
        unsafe
        {
            crocheting::crochet_consts::crochet_consts::set_stitch_dimensions(0.01, 0.01);
        }
        std::fs::create_dir_all("dati_grafici/knitted_test_2").unwrap();
    
        let b =Bezier::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.5, 0.0),
            Some(Point3::new(0.0, 0.0, 0.45)),
            //None,
            None
        );

        b.to_file(Some("dati_grafici/knitted_test_2/curve.mat".to_string()));
        
        
        let s = Cylinder::new(Box::new(b), 0.05, 0.02);
    

        s.to_file(Some("dati_grafici/knitted_test_2/spirale_target.mat".to_string()));

        // Coordinate anelli sovrapposti
        const SAMPLE_ANELLO: i32 = 30;
        let mut my_string = "".to_string();
        for i in 0..s.tau_max_value().floor() as i32
        {
            for n in 0..=30
            {
                let p = s.to_3d_point(
                    TauCoords{tau: i as f64, v: n as f64 / SAMPLE_ANELLO as f64}
                ).unwrap();
                my_string.push_str(&format!("{} {} {}\n", &p.x, &p.y, &p.z));
            }
        }
        
        let mut output = File::create("dati_grafici/knitted_test_2/anelli_sovrapposti.mat").unwrap();
        write!(output, "{}", my_string).unwrap();

        

    
        let mut gr = nodes::CrochetGraph2::new(Box::new(s));
        gr.generate_graph();
        
        gr.start_rows_nodes(Some("dati_grafici/knitted_test_2/startrows.mat".to_string()));

        gr.dot_graph(Some("dati_grafici/knitted_test_2/graph.dot".to_string()));

        gr.to_string_vertices(Some("dati_grafici/knitted_test_2/vertici.mat".to_string()));

        gr.debug_shortrow_level_nodes(Some("dati_grafici/knitted_test_2/vertici_shortrow.mat".to_string()), 1);
        gr.debug_shortrow_level_nodes(Some("dati_grafici/knitted_test_2/vertici_shortrow2.mat".to_string()), 2);

        gr.edge_error_distribution(Some("dati_grafici/knitted_test_2/error_distribution.mat".into()), None);
        gr.metrics_skew(Some("dati_grafici/knitted_test_2/skew_distribution.mat".into()));

        // let mut istr = InstructionsGenerator::new(Box::new(gr));
        // istr.human_instructions(Some("dati_grafici/knitted_test_2/istruzioni.txt".to_string()));
        
    }

    /// Cilindro semplice, dritto
    fn knitted_cylinder()
    {
        unsafe
        {
            crocheting::crochet_consts::crochet_consts::set_stitch_dimensions(0.07, 0.07);
        }
        std::fs::create_dir_all("dati_grafici/knitted_cylinder").unwrap();

        let b =Bezier::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.5, 0.0),
            //Some(Point3::new(0.0, 0.0, 0.45)),
            None,
            None
        );

        b.to_file(Some("dati_grafici/knitted_cylinder/curve.mat".to_string()));
        
        
        let s = Cylinder::new(Box::new(b), 0.11, 0.31);


        s.to_file(Some("dati_grafici/knitted_cylinder/spirale_target.mat".to_string()));


        let mut gr = nodes::CrochetGraph2::new(Box::new(s));
        gr.generate_graph();
        

        gr.dot_graph(Some("dati_grafici/knitted_cylinder/graph.dot".to_string()));

        gr.to_string_vertices(Some("dati_grafici/knitted_cylinder/vertici.mat".to_string()));
        gr.start_rows_nodes(Some("dati_grafici/knitted_cylinder/startrows.mat".to_string()));
        gr.debug_shortrow_level_nodes(Some("dati_grafici/knitted_cylinder/shortrows.mat".into()), 1);

        gr.edge_error_distribution(Some("dati_grafici/knitted_cylinder/error_distribution.mat".into()), None);
        gr.metrics_skew(Some("dati_grafici/knitted_cylinder/skew_distribution.mat".into()));

        // let mut istr = InstructionsGenerator::new(Box::new(gr));
        // istr.human_instructions(Some("dati_grafici/knitted_cylinder/istruzioni.txt".to_string()));
    }

    fn bezier_examples()
    {
        unsafe
        {
            crocheting::crochet_consts::crochet_consts::set_stitch_dimensions(0.1, 0.1);
        }
        let b_linear =Bezier::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.5, 0.0),
            //Some(Point3::new(0.0, 0.0, 0.45)),
            None,
            None
        );
        let b_quad =Bezier::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.5, 0.0),
            Some(Point3::new(0.0, 0.0, 0.45)),
            //None,
            None
        );
        let b_cubic =Bezier::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.5, 0.0),
            Some(Point3::new(0.0, 0.0, 0.45)),
            Some(Point3::new(0.1, 1.0, 2.0)),
        );

        b_linear.to_file(Some("dati_grafici/bezier_examples_curve_linear.mat".to_string()));
        b_quad.to_file(Some("dati_grafici/bezier_examples_curve_quadratic.mat".to_string()));
        b_cubic.to_file(Some("dati_grafici/bezier_examples_curve_cubic.mat".to_string()));
    }


    fn godot_meshviz()
    {
        std::fs::create_dir_all("dati_grafici/godot_meshviz/").unwrap();
        const SCALE: f64 = 33.33333;
        // unsafe{crocheting::crochet_consts::crochet_consts::set_stitch_dimensions(0.006, 0.006)};
        unsafe{crocheting::crochet_consts::crochet_consts::set_stitch_dimensions(0.2, 0.2)};
        let b =Bezier::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 0.1)*SCALE,
            Some(Point3::new(0.0, -0.02, 0.05)*SCALE),
            None
        );
        let s = Cylinder::new(Box::new(b.clone()), 0.05*SCALE, 0.02*SCALE);
        s.to_file(Some("dati_grafici/godot_meshviz/spirale_target.mat".to_string()));
        //let mut gr = CapunamanPaperGraph::new(Box::new(s));
        let mut gr: Box<dyn KnittableGraph> = Box::new(nodes::CapunamanPaperGraph::new(Box::new(s)));
        gr.generate_graph();


        gr.dot_graph(Some("dati_grafici/godot_meshviz/grafogodot.dot".into()));
        gr.to_string_vertices(Some("dati_grafici/godot_meshviz/vertici_capunaman.mat".to_string()));
        gr.edge_error_distribution(Some("dati_grafici/godot_meshviz/error_distribution_capunaman.mat".into()), None);
        gr.edge_error_distribution(
            Some("dati_grafici/godot_meshviz/error_distribution_capunaman_WALE.mat".into()),
            Some(EdgeType::Wale)
        );
        gr.metrics_skew(Some("dati_grafici/godot_meshviz/skew_distribution_capunaman.mat".into()));
        gr.debug_distance_between_rows("dati_grafici/godot_meshviz/distance_between_rows_capunaman.txt");
        let mut st = StitchMesh::new_from_graph(&gr);
        st.array_verts(false);
        st.data_for_blender_viz("dati_grafici/godot_meshviz/blender_facedata_capunaman.txt");
        let istr = InstructionsGenerator::new(&st);
        istr.text_instructions(Some("dati_grafici/godot_meshviz/istruzioni_capunaman.txt"));

        //
        println!("meshviz mio");
        let s = Cylinder::new(Box::new(b), 0.05*SCALE, 0.02*SCALE);
        let mut gr: Box<dyn KnittableGraph> = Box::new(nodes::CrochetGraph2::new(Box::new(s)));
        gr.generate_graph();

        gr.dot_graph(Some("dati_grafici/godot_meshviz/grafogodot.dot".into()));
        gr.to_string_vertices(Some("dati_grafici/godot_meshviz/vertici_gismondi.mat".to_string()));
        gr.debug_shortrow_level_nodes(Some("dati_grafici/godot_meshviz/vertici_gismondi_shortrow1.mat".to_string()), 1);
        gr.debug_shortrow_level_nodes(Some("dati_grafici/godot_meshviz/vertici_gismondi_shortrow2.mat".to_string()), 2);
        gr.edge_error_distribution(Some("dati_grafici/godot_meshviz/error_distribution_gismondi.mat".into()), None);
        gr.edge_error_distribution(
            Some("dati_grafici/godot_meshviz/error_distribution_gismondi_WALE.mat".into()),
            Some(EdgeType::Wale)
        );
        gr.debug_distance_between_rows("dati_grafici/godot_meshviz/distance_between_rows.txt");
        gr.metrics_skew(Some("dati_grafici/godot_meshviz/skew_distribution.mat".into()));

        let mut st = StitchMesh::new_from_graph(&gr);
        st.array_verts(false);
        st.data_for_blender_viz("dati_grafici/godot_meshviz/blender_facedata.txt");
        st.all_face_data_to_file("dati_grafici/godot_meshviz/facedata.txt");

        

        let istr = InstructionsGenerator::new(&st);
        istr.text_instructions(Some("dati_grafici/godot_meshviz/istruzioni.txt"));
    }


    fn capunaman()
    {
        std::fs::create_dir_all("dati_grafici/capunaman/").unwrap();
        unsafe
        {
            crocheting::crochet_consts::crochet_consts::set_stitch_dimensions(0.07, 0.07);
        }

        let b =Bezier::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.5, 0.0),
            //Some(Point3::new(0.0, 0.0, 0.45)),
            None,
            None
        );

        b.to_file(Some("dati_grafici/capunaman/curve.mat".to_string()));
        
        
        let s = Cylinder::new(Box::new(b), 0.11, 0.31);


        s.to_file(Some("dati_grafici/capunaman/spirale_target.mat".to_string()));


        let mut gr = nodes::CapunamanPaperGraph::new(Box::new(s));
        gr.generate_graph();
        

        gr.dot_graph(Some("dati_grafici/capunaman/graph.dot".to_string()));

        gr.to_string_vertices(Some("dati_grafici/capunaman/vertici.mat".to_string()));
        gr.start_rows_nodes(Some("dati_grafici/capunaman/startrows.mat".to_string()));

        gr.edge_error_distribution(Some("dati_grafici/capunaman/error_distribution.mat".into()), None);
        gr.metrics_skew(Some("dati_grafici/capunaman/skew_distribution.mat".into()));

        // let mut istr = InstructionsGenerator::new(Box::new(gr));
        // istr.human_instructions(Some("dati_grafici/capunaman/istruzioni.txt".to_string()));
    }


    fn capunaman_knitted_test_2()
    {
        std::fs::create_dir_all("dati_grafici/capunaman_knitted_test_2/").unwrap();
        unsafe
        {
            crocheting::crochet_consts::crochet_consts::set_stitch_dimensions(0.03, 0.03);
        }

        let b =Bezier::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.5, 0.0),
            Some(Point3::new(0.0, 0.0, 0.45)),
            //None,
            None
        );

        b.to_file(Some("dati_grafici/capunaman_knitted_test_2/curve.mat".to_string()));
        
        
        let s = Cylinder::new(Box::new(b), 0.05, 0.02);
    

        s.to_file(Some("dati_grafici/capunaman_knitted_test_2/spirale_target.mat".to_string()));

    
        let mut gr = nodes::CapunamanPaperGraph::new(Box::new(s));
        gr.generate_graph();
        

        gr.dot_graph(Some("dati_grafici/capunaman_knitted_test_2/_graph.dot".to_string()));

        gr.to_string_vertices(Some("dati_grafici/capunaman_knitted_test_2/vertici.mat".to_string()));

        gr.edge_error_distribution(Some("dati_grafici/capunaman_knitted_test_2/error_distribution.mat".into()), None);
        gr.metrics_skew(Some("dati_grafici/capunaman_knitted_test_2/skew_distribution.mat".into()));

        // let mut istr = InstructionsGenerator::new(Box::new(gr));
        // istr.human_instructions(Some("dati_grafici/capunaman_knitted_test_2/istruzioni.txt".to_string()));
    }

    fn multishape_1_crochetgraph()
    {
        std::fs::create_dir_all("dati_grafici/multishape_1_crochetgraph/").unwrap();
        unsafe
        {
            crocheting::crochet_consts::crochet_consts::set_stitch_dimensions(0.005, 0.005);
        }

        let b =Bezier::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 0.2),
            None,
            None
        );

        b.to_file(Some("dati_grafici/multishape_1_crochetgraph/curve.mat".to_string()));
        
        let base = Circle::new(0.02);
        let first = MultiPointShape::new(
            vec![
                (0.0, Point3::new(0.07, 0.0, 0.0)),
                (0.25, Point3::new(0.07, 0.07, 0.0)),
                (0.4, Point3::new(0.0, 0.07, 0.0)),
                (0.5, Point3::new(-0.02, 0.0, 0.0)),
                (0.75, Point3::new(0.0, -0.02, 0.0)),
            ]
        );
        let last = MultiPointShape::new(
            vec![
                (0.0, Point3::new(0.02, 0.0, 0.0)),
                (0.25, Point3::new(0.0, 0.02, 0.0)),
                (0.5, Point3::new(-0.02, 0.0, 0.0)),
                (0.75, Point3::new(0.0, -0.02, 0.0)),
            ]
        );
        
        let s = MultiShapeSurface::new(
            Box::new(b.clone()),
            vec![
                (0.0, Box::new(base)),
                (0.3, Box::new(first.clone())),
                (1.0, Box::new(last.clone()))
            ]
        );


        s.to_file(Some("dati_grafici/multishape_1_crochetgraph/spirale_target.mat".to_string()));



        unsafe {
            crocheting::crochet_consts::crochet_consts::set_flag_generate_short_row_flag(false);
        }
        let mut gr: Box<dyn KnittableGraph> = Box::new(nodes::CrochetGraph2::new(Box::new(s)));
        gr.generate_graph();
        gr.to_string_vertices(Some("dati_grafici/multishape_1_crochetgraph/vertici_no_shortrow_generation.mat".to_string()));

        gr.debug_distance_between_rows("dati_grafici/multishape_1_crochetgraph/distance_between_rows_NOGENERATE.txt");

        unsafe {
            crocheting::crochet_consts::crochet_consts::set_flag_generate_short_row_flag(true);
        }
        let s = MultiShapeSurface::new(
            Box::new(b),
            vec![
                (0.0, Box::new(base)),
                (0.3, Box::new(first)),
                (1.0, Box::new(last))
            ]
        );
        let mut gr: Box<dyn KnittableGraph> = Box::new(nodes::CrochetGraph2::new(Box::new(s)));
        gr.generate_graph();
        gr.to_string_vertices(Some("dati_grafici/multishape_1_crochetgraph/vertici.mat".to_string()));

        gr.dot_graph(Some("dati_grafici/multishape_1_crochetgraph/graph.dot".to_string()));

        gr.start_rows_nodes(Some("dati_grafici/multishape_1_crochetgraph/startrows.mat".to_string()));
        gr.debug_shortrow_level_nodes(Some("dati_grafici/multishape_1_crochetgraph/shortrows.mat".into()), 1);
        gr.debug_shortrow_level_nodes(Some("dati_grafici/multishape_1_crochetgraph/shortrows2.mat".into()), 2);

        gr.edge_error_distribution(Some("dati_grafici/multishape_1_crochetgraph/error_distribution.mat".into()), None);
        gr.metrics_skew(Some("dati_grafici/multishape_1_crochetgraph/skew_distribution.mat".into()));
        gr.metrics_knitting_direction(Some("dati_grafici/multishape_1_crochetgraph/metrics_knitting_dir.mat".into()));
        gr.debug_distance_between_rows("dati_grafici/multishape_1_crochetgraph/distance_between_rows.txt");

        let st = StitchMesh::new_from_graph(&gr);

        let istr = InstructionsGenerator::new(&st);
        istr.text_instructions(Some("dati_grafici/multishape_1_crochetgraph/istruzioni.txt"));
    }


    fn multishape_1_capunaman()
    {
        std::fs::create_dir_all("dati_grafici/multishape_1_capunaman/").unwrap();
        unsafe
        {
            crocheting::crochet_consts::crochet_consts::set_stitch_dimensions(0.005, 0.005);
        }

        let b =Bezier::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 0.2),
            None,
            None
        );

        b.to_file(Some("dati_grafici/multishape_1_capunaman/curve.mat".to_string()));
        
        let base = Circle::new(0.02);
        let first = MultiPointShape::new(
            vec![
                (0.0, Point3::new(0.07, 0.0, 0.0)),
                (0.25, Point3::new(0.07, 0.07, 0.0)),
                (0.4, Point3::new(0.0, 0.07, 0.0)),
                (0.5, Point3::new(-0.02, 0.0, 0.0)),
                (0.75, Point3::new(0.0, -0.02, 0.0)),
            ]
        );
        let last = MultiPointShape::new(
            vec![
                (0.0, Point3::new(0.02, 0.0, 0.0)),
                (0.25, Point3::new(0.0, 0.02, 0.0)),
                (0.5, Point3::new(-0.02, 0.0, 0.0)),
                (0.75, Point3::new(0.0, -0.02, 0.0)),
            ]
        );
        
        let s = MultiShapeSurface::new(
            Box::new(b),
            vec![
                (0.0, Box::new(base)),
                (0.3, Box::new(first)),
                (1.0, Box::new(last))
            ]
        );


        s.to_file(Some("dati_grafici/multishape_1_capunaman/spirale_target.mat".to_string()));


        let mut gr = nodes::CapunamanPaperGraph::new(Box::new(s));
        gr.generate_graph();
        

        gr.dot_graph(Some("dati_grafici/multishape_1_capunaman/graph.dot".to_string()));

        gr.to_string_vertices(Some("dati_grafici/multishape_1_capunaman/vertici.mat".to_string()));
        gr.start_rows_nodes(Some("dati_grafici/multishape_1_capunaman/startrows.mat".to_string()));
        gr.debug_shortrow_level_nodes(Some("dati_grafici/multishape_1_capunaman/shortrows.mat".into()), 1);

        gr.edge_error_distribution(Some("dati_grafici/multishape_1_capunaman/error_distribution.mat".into()), None);
        gr.metrics_skew(Some("dati_grafici/multishape_1_capunaman/skew_distribution.mat".into()));
        gr.metrics_knitting_direction(Some("dati_grafici/multishape_1_capunaman/metrics_knitting_dir.mat".into()));
        gr.debug_distance_between_rows("dati_grafici/multishape_1_capunaman/distance_between_rows.txt");

        // let mut istr = InstructionsGenerator::new(Box::new(gr));
        // istr.human_instructions(Some("dati_grafici/multishape_1_capunaman/istruzioni.txt".to_string()));
    }


    fn oca_celeste()
    {
        std::fs::create_dir_all("dati_grafici/oca_celeste/").unwrap();
        unsafe
        {
            crocheting::crochet_consts::crochet_consts::set_stitch_dimensions(0.45, 0.45);
        }
        let b =Bezier::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 14.0),
            None,
            None
        );

        let forme = vec![
            (0.0,   MultiPointShape::new_ellipse(0.7, 0.7, (0.0, 0.0))),
            (2.5,   MultiPointShape::new_ellipse(4.0, 6.0, (0.0, 0.0))),
            (4.0,   MultiPointShape::new_ellipse(5.0, 7.5, (0.0, 0.0))),
            (5.0,   MultiPointShape::new_ellipse(5.0, 7.5, (0.0, 0.0))),
            (7.0,   MultiPointShape::new_ellipse(2.5, 2.5, (0.0, 4.0))),
            (8.5,   MultiPointShape::new_ellipse(1.0, 1.0, (0.0, 5.0))),
            (10.0,  MultiPointShape::new_ellipse(1.0, 1.0, (0.0, 4.75))),
            (12.0,  MultiPointShape::new_ellipse(1.0, 1.0, (0.0, 4.5))),
            (13.0,  MultiPointShape::new_ellipse(1.5, 1.5, (0.0, 4.5))),
            (14.0,  MultiPointShape::new_ellipse(1.0, 1.0, (0.0, 4.5))),
        ];

        let mut forme_2: Vec<(f64, Box<dyn Shape>)> = Vec::new();
        for f in forme
        {
            forme_2.push(
                (f.0 / 14.0, Box::new(f.1))
            );
        }

        let s = MultiShapeSurface::new(Box::new(b), forme_2);
        s.to_file(Some("dati_grafici/oca_celeste/spirale_target.mat".to_string()));

        unsafe
        {
            crocheting::crochet_consts::crochet_consts::set_flag_generate_short_row_flag(true);
        }
        let mut gr: Box<dyn KnittableGraph> = Box::new(nodes::CrochetGraph2::new(Box::new(s.clone())));
        gr.generate_graph();
        gr.debug_distance_between_rows("dati_grafici/oca_celeste/distance_between_rows.txt");
        gr.to_string_vertices(Some("dati_grafici/oca_celeste/vertici.mat".to_string()));
        gr.start_rows_nodes(Some("dati_grafici/oca_celeste/startrows.mat".to_string()));
        gr.debug_shortrow_level_nodes(Some("dati_grafici/oca_celeste/shortrows.mat".into()), 1);
        gr.debug_shortrow_level_nodes(Some("dati_grafici/oca_celeste/shortrows2.mat".into()), 2);
        gr.debug_shortrow_level_nodes(Some("dati_grafici/oca_celeste/shortrows3.mat".into()), 3);

        gr.dot_graph(Some("dati_grafici/oca_celeste/graph.dot".to_string()));

        

        gr.edge_error_distribution(Some("dati_grafici/oca_celeste/error_distribution.mat".into()), None);
        gr.edge_error_distribution(Some("dati_grafici/oca_celeste/error_distribution_wale.mat".into()), Some(EdgeType::Wale));
        gr.metrics_skew(Some("dati_grafici/oca_celeste/skew_distribution.mat".into()));
        gr.metrics_knitting_direction(Some("dati_grafici/oca_celeste/metrics_knitting_dir.mat".into()));

        let mut st = StitchMesh::new_from_graph(&gr);
        st.array_verts(false);
        st.data_for_blender_viz("dati_grafici/oca_celeste/blender_data.txt");
        st.all_face_data_to_file("dati_grafici/oca_celeste/facedata.txt");

        let istr = InstructionsGenerator::new(&st);
        istr.text_instructions(Some("dati_grafici/oca_celeste/istruzioni.txt"));



        // Versione capunaman
        let mut gr: Box<dyn KnittableGraph> = Box::new(nodes::CapunamanPaperGraph::new(Box::new(s)));
        gr.generate_graph();
        gr.dot_graph(Some("dati_grafici/oca_celeste/graph_capunaman.dot".to_string()));
        gr.to_string_vertices(Some("dati_grafici/oca_celeste/vertici_capunaman.mat".to_string()));
        gr.edge_error_distribution(Some("dati_grafici/oca_celeste/error_distribution_capunaman.mat".into()), None);
        gr.edge_error_distribution(Some("dati_grafici/oca_celeste/error_distribution_wale_capunaman.mat".into()), Some(EdgeType::Wale));
        gr.metrics_skew(Some("dati_grafici/oca_celeste/skew_distribution_capunaman.mat".into()));
        gr.metrics_knitting_direction(Some("dati_grafici/oca_celeste/metrics_knitting_dir_capunaman.mat".into()));

        let mut st = StitchMesh::new_from_graph(&gr);
        st.array_verts(false);
        st.data_for_blender_viz("dati_grafici/oca_celeste/capunaman_blender_data.txt");

        unsafe
        {
            crocheting::crochet_consts::crochet_consts::set_flag_generate_short_row_flag(true);
        }
    }

    fn beastie()
    {
        std::fs::create_dir_all("dati_grafici/beastie/").unwrap();
        unsafe
        {
            crocheting::crochet_consts::crochet_consts::set_stitch_dimensions(0.8, 0.8);// Centimetri
        }
        let s = modelli::modelli::beastie();
        s.to_file(Some("dati_grafici/beastie/spirale_target.mat".to_string()));

        let mut gr: Box<dyn KnittableGraph> = Box::new(nodes::CrochetGraph2::new(s));
        gr.generate_graph();
        gr.dot_graph(Some("dati_grafici/beastie/graph.dot".to_string()));
        gr.to_string_vertices(Some("dati_grafici/beastie/vertici.mat".to_string()));
        gr.start_rows_nodes(Some("dati_grafici/beastie/startrows.mat".to_string()));
        gr.debug_shortrow_level_nodes(Some("dati_grafici/beastie/shortrows.mat".into()), 1);

        //gr.dot_graph(Some("tmp_graph.dot".into()));

        gr.edge_error_distribution(Some("dati_grafici/beastie/error_distribution.mat".into()), None);
        gr.metrics_skew(Some("dati_grafici/beastie/skew_distribution.mat".into()));

        let st = StitchMesh::new_from_graph(&gr);
        let istr = InstructionsGenerator::new(&st);
        istr.text_instructions(Some("dati_grafici/beastie/istruzioni.txt"));
    }

    fn spiral_check()
    {
        std::fs::create_dir_all("dati_grafici/spiral_check/").unwrap();
        unsafe
        {
            crocheting::crochet_consts::crochet_consts::set_stitch_dimensions(0.5, 0.5);// Centimetri
        }
        let b =Bezier::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 8.0),
            None,
            None
        );
        let f1 = MultiPointShape::new(
            vec![
                (0.0,  Point3::new(1.0, 0.0, -2.0)),
                (0.125, Point3::new(1.0, 1.0, -1.0)),
                (0.25,  Point3::new(0.0, 1.0, 0.0)),
                (0.325, Point3::new(-1.0, 1.0, -1.0)),
                (0.5,   Point3::new(-1.0, 0.0, -2.0)),
                (0.625, Point3::new(-1.0, -1.0, -1.0)),
                (0.75,  Point3::new(0.0, -1.0, 0.0)),
                (0.825, Point3::new(1.0, -1.0, -1.0)),
            ]
        );
        let f2 = MultiPointShape::new_magic_ring((0.0, 0.0));

        let s = MultiShapeSurface::new(
            Box::new(b),
            vec![
                (0.0, Box::new(f1)),
                (1.0, Box::new(f2))
            ]
        );
        
        let subdivision_tau: usize = 25;
        const SUBDIVISION_V: usize = 50;
        let mut file_output = File::create("dati_grafici/spiral_check/superficie.mat").unwrap();
        // Genero prima i vertici
        for i in 0..=subdivision_tau
        {
            let tau = (i as f64 / subdivision_tau as f64)*s.tau_max_value();
            for j in 0..SUBDIVISION_V
            {
                let v = j as f64 / SUBDIVISION_V as f64;

                let p = s.to_3d_point(TauCoords{tau: tau, v: v}).unwrap();
                write!(file_output, "{} {} {}\n", p.x, p.y, p.z);
            }
        }

        let mut gr_1 = nodes::CapunamanPaperGraph::new(Box::new(s.clone()));
        gr_1.generate_graph();
        gr_1.to_string_vertices(Some("dati_grafici/spiral_check/spirale.mat".into()));
    }


    fn plots_for_compliance()
    {
        std::fs::create_dir_all("dati_grafici/plots_for_compliance/").unwrap();
        unsafe
        {
            crocheting::crochet_consts::crochet_consts::set_stitch_dimensions(0.4, 0.4);// Centimetri
        }

        let b =Bezier::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 5.0),
            None,
            None
        );

        let s = Cylinder::new(Box::new(b), 3.0, 1.0);

        let mut gr_1 = nodes::CapunamanPaperGraph::new(Box::new(s.clone()));
        gr_1.generate_graph();
        gr_1.to_string_vertices(Some("dati_grafici/plots_for_compliance/ogni_riga.mat".into()));
        gr_1.start_rows_nodes(Some("dati_grafici/plots_for_compliance/startrows.mat".into()));

        let mut gr_2 = nodes::CapunamanPaperGraph::new(Box::new(s));
        gr_2.generate_graph();
        gr_2.to_string_vertices(Some("dati_grafici/plots_for_compliance/curva_continua.mat".into()));
        gr_2.start_rows_nodes(Some("dati_grafici/plots_for_compliance/startrows_continua.mat".into()));
    }



    fn time_benchmark_godot_meshviz()
    {
        println!("godotmeshviz");
        let b =Bezier::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 0.1),
            Some(Point3::new(0.0, -0.02, 0.05)),
            None
        );

        

        for (i, dim) in [0.02, 0.018, 0.015, 0.01, 0.005, 0.004, 0.0035, 0.003, 0.0028, 0.0025, 0.001].iter().enumerate()
        {
            unsafe
            {
                crocheting::crochet_consts::crochet_consts::set_stitch_dimensions(*dim, *dim);// Centimetri
            }
            let s: Box<dyn UVSurface> = Box::new(Cylinder::new(Box::new(b.clone()), 0.05, 0.02));

            static_run_benchmark(s, "godot_meshviz", i != 0);
        }
    }

    fn time_benchmark_oca_celeste()
    {
        println!("oca_celeste");
        let b =Bezier::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 14.0),
            None,
            None
        );

        let forme = vec![
            (0.0,   MultiPointShape::new_ellipse(0.7, 0.7, (0.0, 0.0))),
            (2.5,   MultiPointShape::new_ellipse(4.0, 6.0, (0.0, 0.0))),
            (4.0,   MultiPointShape::new_ellipse(5.0, 7.5, (0.0, 0.0))),
            (5.0,   MultiPointShape::new_ellipse(5.0, 7.5, (0.0, 0.0))),
            (7.0,   MultiPointShape::new_ellipse(2.5, 2.5, (0.0, 4.0))),
            (8.5,   MultiPointShape::new_ellipse(1.0, 1.0, (0.0, 5.0))),
            (10.0,  MultiPointShape::new_ellipse(1.0, 1.0, (0.0, 4.75))),
            (12.0,  MultiPointShape::new_ellipse(1.0, 1.0, (0.0, 4.5))),
            (13.0,  MultiPointShape::new_ellipse(1.5, 1.5, (0.0, 4.5))),
            (14.0,  MultiPointShape::new_ellipse(1.0, 1.0, (0.0, 4.5))),
        ];

        let mut forme_2: Vec<(f64, Box<dyn Shape>)> = Vec::new();
        for f in forme
        {
            forme_2.push(
                (f.0 / 14.0, Box::new(f.1))
            );
        }

        for (i, dim) in [0.5, 0.4, 0.3, 0.25, 0.2, 0.16, 0.125, 0.1].iter().enumerate()
        {
            unsafe
            {
                crocheting::crochet_consts::crochet_consts::set_stitch_dimensions(*dim, *dim);// Centimetri
            }
            let s: Box<dyn UVSurface> = Box::new(MultiShapeSurface::new(Box::new(b.clone()), forme_2.clone()));

            static_run_benchmark(s, "oca_celeste", i != 0);
        }
    }


    // Calcola per le quattro superfici del paper
    // lo stiramento e skew per ogni risoluzione data
    fn mean_variance_over_stitch_size()
    {
        const SCALE_C: f64 = 33.3333333;
        const SCALE_D: f64 = 16.0*1.25;
        const DIM_STITCH_C: [f64; 3] = [0.25, 0.2, 0.15];
        const DIM_STITCH_D: [f64; 3] = [0.25, 0.2, 0.15];

        let curve_c =Bezier::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 0.1)*SCALE_C,
            Some(Point3::new(0.0, -0.02, 0.05)*SCALE_C),
            None
        );

        let curve_d =Bezier::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.5, 0.0)*SCALE_D,
            Some(Point3::new(0.0, 0.0, 0.45)*SCALE_D),
            //None,
            None
        );

        // Per ogni dim stitch, genero i dati di gismondi e capunaman
        for (i, dim) in DIM_STITCH_C.iter().enumerate()
        {
            unsafe
            {
                crocheting::crochet_consts::crochet_consts::set_stitch_dimensions(*dim, *dim);// Centimetri
            }
            println!("Sup_c a {}", dim);
            let sup_c: Box<dyn UVSurface> = Box::new(Cylinder::new(Box::new(curve_c.clone()), 0.05*SCALE_C, 0.02*SCALE_C));
            static_run_benchmark(sup_c, "sup_c", i != 0);
        }


        for (i, dim) in DIM_STITCH_D.iter().enumerate()
        {
            unsafe
            {
                crocheting::crochet_consts::crochet_consts::set_stitch_dimensions(*dim, *dim);// Centimetri
            }
            println!("Sup_d a {}", dim);
            let sup_d: Box<dyn UVSurface> = Box::new(Cylinder::new(Box::new(curve_d.clone()), 0.05*SCALE_D, 0.02*SCALE_D));
            static_run_benchmark(sup_d, "sup_d", i != 0);
        }
    }


    fn static_run_benchmark(surf: Box<dyn UVSurface>, nome: &str, append: bool)
    {
        std::fs::create_dir_all("dati_grafici/time_benchmark/").unwrap();
        std::fs::create_dir_all("dati_grafici/paper_results_mean_variance/").unwrap();
        std::fs::create_dir_all("dati_grafici/paper_results_mean_variance/sup_c/").unwrap();
        std::fs::create_dir_all("dati_grafici/paper_results_mean_variance/sup_d/").unwrap();

        let current_dim = crocheting::crochet_consts::crochet_consts::get_width();

        let mut file_capunaman: std::fs::File;
        let mut file_gismondi: std::fs::File;
        if append
        {
            file_capunaman = fs::OpenOptions::new()
            .write(true)
            .append(append)
            .open(
                format!("dati_grafici/time_benchmark/{}_capunaman.mat", nome)
            )
            .unwrap();
            file_gismondi = fs::OpenOptions::new()
            .write(true)
            .append(append)
            .open(
                format!("dati_grafici/time_benchmark/{}_gismondi.mat", nome)
            )
            .unwrap();
        }
        else
        {
            file_capunaman = File::create(format!("dati_grafici/time_benchmark/{}_capunaman.mat", nome)).unwrap();
            file_gismondi = File::create(format!("dati_grafici/time_benchmark/{}_gismondi.mat", nome)).unwrap();
        }
        {
            println!("\nlanciato a dimensione {}", crocheting::crochet_consts::crochet_consts::get_height());
            println!("   Capunaman");
            {
                let start = Instant::now();
                let mut gr: Box<dyn KnittableGraph> = Box::new(nodes::CapunamanPaperGraph::new(surf.clone()));
                println!("------inizial");
                gr.generate_graph();
                println!("------Generato");
                let fine_calcolo_reticolo = Instant::now();
                
                let avvio_render_vertici = Instant::now();
                
                let end = Instant::now();
                gr.edge_error_distribution(
                    Some(
                        format!("dati_grafici/paper_results_mean_variance/{}/{}_capunaman_edge_error.txt", nome, current_dim)
                    ),
                    None);
                gr.metrics_skew(
                    Some(
                        format!("dati_grafici/paper_results_mean_variance/{}/{}_capunaman_skew.txt", nome, current_dim)
                    )
                );

                let mut st = StitchMesh::new_from_graph(&gr);
                let istr = InstructionsGenerator::new(&st);
                istr.text_instructions(Some(
                    &format!("../cannavo_dati_integrativi/istruzioni/{}_{}_capunaman.txt", nome, current_dim)
                ));

                let dati_capunaman: [Duration; 3] = [
                    fine_calcolo_reticolo.duration_since(start),
                    avvio_render_vertici.duration_since(fine_calcolo_reticolo),
                    end.duration_since(avvio_render_vertici)
                ];


                let mut time_file = File::create(format!("dati_grafici/paper_results_mean_variance/{}/{}_time_capunaman.txt", nome, current_dim)).unwrap();
                writeln!(time_file, "{}", dati_capunaman[0].as_millis());
                {
                    write!(file_capunaman, "{} {} {} {}\n",
                        gr.verts_count(),
                        dati_capunaman[0].as_millis(),
                        dati_capunaman[1].as_millis(),
                        dati_capunaman[2].as_millis()).unwrap();
                }
            }
            
            println!("   Gismondi");
            {
                let start = Instant::now();
                let mut gr: Box<dyn KnittableGraph> = Box::new(nodes::CrochetGraph2::new(surf.clone()));
                gr.generate_graph();
                let fine_calcolo_reticolo = Instant::now();
                // let mut st = StitchMesh::new_from_graph(&gr);
                
                let avvio_render_vertici = Instant::now();
                //st.array_verts();
                let end = Instant::now();
    
                gr.edge_error_distribution(
                    Some(
                        format!("dati_grafici/paper_results_mean_variance/{}/{}_gismondi_edge_error.txt", nome, current_dim)
                    ),
                    None);
                gr.metrics_skew(
                    Some(
                        format!("dati_grafici/paper_results_mean_variance/{}/{}_gismondi_skew.txt",nome,  current_dim)
                    )
                );


                let mut st = StitchMesh::new_from_graph(&gr);
                let istr = InstructionsGenerator::new(&st);
                istr.text_instructions(Some(
                    &format!("../cannavo_dati_integrativi/istruzioni/{}_{}_gismondi.txt", nome, current_dim)
                ));

                let dati_gismondi: [Duration; 3] = [
                    fine_calcolo_reticolo.duration_since(start),
                    avvio_render_vertici.duration_since(fine_calcolo_reticolo),
                    end.duration_since(avvio_render_vertici)
                ];
                let mut time_file = File::create(format!("dati_grafici/paper_results_mean_variance/{}/{}_time_gismondi.txt", nome, current_dim)).unwrap();
                writeln!(time_file, "{}", dati_gismondi[0].as_millis());
                {
                    write!(file_gismondi, "{} {} {} {}\n",
                        gr.verts_count(),
                        dati_gismondi[0].as_millis(),
                        dati_gismondi[1].as_millis(),
                        dati_gismondi[2].as_millis()).unwrap();
                }
            }
            
        }
    }



    fn asymmetric_expansion()
    {
        std::fs::create_dir_all("dati_grafici/asymmetric_expansion/").unwrap();
        unsafe
        {
            crocheting::crochet_consts::crochet_consts::set_stitch_dimensions(0.4, 0.4);// Centimetri
        }

        let b =Bezier::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 5.0),
            None,
            None
        );

        let s = MultiShapeSurface::new(Box::new(b), 
            vec![
                (0.0, Box::new(MultiPointShape::new_magic_ring((0.0, 0.0)))),
                (0.5, Box::new(
                    MultiPointShape::new_ellipse(3.0, 8.0, (0.0, 0.0))
                )),
                (1.0, Box::new(
                    MultiPointShape::new_ellipse(3.0, 8.0, (0.0, 0.0))
                ))
            ]
        );

        let mut gr_1 = nodes::CrochetGraph2::new(Box::new(s.clone()));
        gr_1.generate_graph();
        gr_1.to_string_vertices(Some("dati_grafici/asymmetric_expansion/gismondi:vertici.mat".into()));
        gr_1.edge_error_distribution(Some("dati_grafici/asymmetric_expansion/error_distribution.mat".into()), None);
        gr_1.metrics_skew(Some("dati_grafici/asymmetric_expansion/skew_distribution.mat".into()));


        let mut gr_2 = nodes::CapunamanPaperGraph::new(Box::new(s));
        gr_2.generate_graph();
        gr_2.to_string_vertices(Some("dati_grafici/asymmetric_expansion/curva_continua.mat".into()));
        gr_2.start_rows_nodes(Some("dati_grafici/asymmetric_expansion/startrows_continua.mat".into()));
    }
}