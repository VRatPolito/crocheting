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

pub mod nodes {
    use crate::crochet_consts::crochet_consts;
    use crate::errors::errors::CrochetError;
    use crate::param_surface::surface::TauCoords;
    use crate::param_surface::surface::UVSurface;
    use nalgebra::{Point3, Vector3};
    use petgraph::dot::Dot;
    use petgraph::graph::{EdgeIndex, Graph, NodeIndex};
    use petgraph::visit::EdgeRef;
    use petgraph::visit::IntoNodeReferences;
    use petgraph::Direction;
    use std::collections::HashSet;
    use std::f64::consts::PI;
    use std::fs::File;
    use std::io::Write;

    #[derive(Debug, PartialEq, Clone)]
    pub enum EdgeType {
        Wale,                    //Colonna verticale
        Course,                  //Laterale
        CourseNextOverrideStart, // Salto, inizio shortrow
        CourseNextOverrideEnd,   //Quando devo saltare, fine shortrow
    }
    impl std::fmt::Display for EdgeType {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            write!(f, "{}", format!("{:?}", self))
        }
    }

    #[derive(Clone)]
    pub struct CrochetVertex {
        coords: TauCoords,
        short_row_level: i32,
        target_tau: Option<f64>, //position_3d: Point3<f64>//una sorta di cache per non ricalcolarla ogni volta. Da aggiungere dopo per ottimizzare
    }
    impl std::fmt::Debug for CrochetVertex {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            write!(
                f,
                "{}",
                format!(
                    "{:.3} | {:.3}, shortrow: {}. TARGET: {:.3}
                    ",
                    self.coords.tau,
                    self.coords.v,
                    self.short_row_level,
                    self.target_tau.unwrap()
                )
            )
        }
    }
    impl CrochetVertex {
        fn height_target(&self) -> TauCoords {
            return TauCoords {
                tau: self.target_tau.unwrap(),
                v: self.target_tau.unwrap() - self.target_tau.unwrap().floor(),
            };
        }

        fn move_right_on_helix(&self, delta: f64) -> CrochetVertex {
            let mut new_coords = self.coords.clone();
            new_coords = new_coords
                + TauCoords {
                    tau: delta,
                    v: delta,
                };

            let mut new_tau: Option<f64>;
            if let Some(e) = self.target_tau {
                new_tau = Some(e + delta);
            } else {
                new_tau = None;
            }

            if new_coords.tau < 0.0 {
                new_coords.tau = 0.0;
                new_tau = Some(1.0);
            }
            #[cfg(debug_assertions)]
            {
                // println!("Move on helix. Sposto {:.3?} di {:.5}. Ottengo {:.3?}", past_coords, delta, new_coords);
            }

            CrochetVertex {
                coords: new_coords,
                short_row_level: self.short_row_level,
                target_tau: new_tau,
            }
        }
    }

    pub struct CrochetEdge {
        pub link_type: EdgeType,
    }
    impl std::fmt::Debug for CrochetEdge {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            write!(f, "{}", format!("{:?}", self.link_type))
        }
    }

    ///
    /// ```
    /// use nalgebra::{Point3};
    /// use crocheting::nodes::nodes::internal_compute_skew;
    /// assert_eq!(
    ///     internal_compute_skew(
    ///     [
    ///         Point3::new(0.0, 0.0, 0.0),
    ///         Point3::new(2.0, 0.0, 0.0),
    ///     ],
    ///     [
    ///         Point3::new(0.0, 3.0, 0.0),
    ///         Point3::new(2.0, 3.0, 0.0),
    ///     ]
    ///     ),
    ///     0.0
    /// );
    ///
    /// // In diagonale
    /// assert_eq!(
    ///     internal_compute_skew(
    ///     [
    ///         Point3::new(0.0, 0.0, 0.0),
    ///         Point3::new(3.0, 3.0, 0.0),
    ///     ],
    ///     [
    ///         Point3::new(1.0, 1.0, 1.0),
    ///         Point3::new(2.0, 2.0, 1.0),
    ///     ]
    ///     ),
    ///     0.0
    /// );
    /// ```
    pub fn internal_compute_skew(down: [Point3<f64>; 2], up: [Point3<f64>; 2]) -> f64 {
        let vertical_dir: Vector3<f64> = (Point3::new(
            (up[0].x + up[1].x) / 2.0,
            (up[0].y + up[1].y) / 2.0,
            (up[0].z + up[1].z) / 2.0,
        ) - Point3::new(
            (down[0].x + down[1].x) / 2.0,
            (down[0].y + down[1].y) / 2.0,
            (down[0].z + down[1].z) / 2.0,
        ))
        .normalize();

        let horizontal_dir: Vector3<f64> = (Point3::new(
            (down[1].x + up[1].x) / 2.0,
            (down[1].y + up[1].y) / 2.0,
            (down[1].z + up[1].z) / 2.0,
        ) - Point3::new(
            (down[0].x + up[0].x) / 2.0,
            (down[0].y + up[0].y) / 2.0,
            (down[0].z + up[0].z) / 2.0,
        ))
        .normalize();

        // let teta_e: f64;
        // if nalgebra::distance_squared(&up[0], &up[1]) < 0.000000001
        // {
        //     teta_e = PI/3.0;
        // }
        // else
        // {
        //     teta_e = PI/2.0;
        // }
        // let teta_max = tetas[0].max(tetas[1]);
        // let teta_min = tetas[0].min(tetas[1]);

        // let

        return vertical_dir.dot(&horizontal_dir).abs();
    }

    pub trait KnittableGraph {
        // Default implementation in traits
        // https://users.rust-lang.org/t/best-practices-when-defining-a-default-implementation-for-a-traits-method/2033/2

        /// Funzione da chiamare dopo avere costruito la superficie
        /// per generare il grafo interno.
        /// Qui risiede l'algoritmo principale
        ///
        /// Chiamare questa prima di provare a generare le istruzioni o la visualizzazione
        fn generate_graph(&mut self);

        /// Accesso al grafo sottostante in sola lettura
        fn get_graph<'a>(&'a self) -> &'a Graph<CrochetVertex, CrochetEdge>;

        /// Accesso alla classe che eredita da UVSurface in sola lettura
        fn get_surface<'a>(&'a self) -> &'a Box<dyn UVSurface>;

        /// Il nodo iniziale, utile per inizializzare l'algoritmo di tracing
        /// o di visualizzazione
        fn get_start_node(&self) -> NodeIndex;

        fn edge_length(&self, edge: petgraph::prelude::EdgeIndex) -> f64 {
            let start_end_coords = self.get_graph().edge_endpoints(edge).unwrap();

            #[cfg(debug_assertions)]
            {
                println!("{}", self.get_surface().tau_max_value());
            }

            return self
                .get_surface()
                .geodesic_dist(
                    self.get_graph()
                        .node_weight(start_end_coords.0)
                        .unwrap()
                        .coords,
                    self.get_graph()
                        .node_weight(start_end_coords.1)
                        .unwrap()
                        .coords,
                    None,
                )
                .unwrap();
        }

        fn shortrows_limits(&self) -> Vec<(NodeIndex, NodeIndex)> {
            let mut found_nodes: Vec<(NodeIndex, NodeIndex)> = Vec::new();
            for ed in self.get_graph().edge_references() {
                if ed.weight().link_type == EdgeType::CourseNextOverrideStart {
                    let original_start_point = self.find_next_course(ed.source(), false).unwrap();
                    // Riga dispari
                    {
                        let start_point = original_start_point;
                        let mut end_point: NodeIndex = start_point;
                        loop {
                            let next: Option<NodeIndex> = self.find_next_course(end_point, false);
                            if next.is_none() {
                                break;
                            }
                            end_point = next.unwrap();
                        }
                        found_nodes.push((start_point, end_point));
                    }

                    // Riga pari
                    {
                        let start_point = self.find_next_wale(original_start_point).unwrap();
                        let mut end_point: NodeIndex = start_point;
                        loop {
                            #[cfg(debug_assertions)]
                            {
                                println!("{}", self.print_node(end_point));
                            }
                            let next: NodeIndex = self.find_next_course(end_point, false).unwrap();
                            end_point = next;

                            if self
                                .get_graph()
                                .edges_directed(next, Direction::Outgoing)
                                .filter(|e| e.weight().link_type == EdgeType::CourseNextOverrideEnd)
                                .count()
                                > 0
                            {
                                break;
                            }
                        }
                        found_nodes.push((start_point, end_point));
                    }
                }
            }
            return found_nodes;
        }

        /// Per ogni nodo, calcola la
        /// distanza rispetto alla spirale target
        fn debug_distance_between_rows(&self, filename: &str) {
            let mut cur_node = self.get_start_node();
            let mut outp = "".to_string();
            while let Some(nuovo) = self.find_next_course(cur_node, true) {
                let w = self.get_graph().node_weight(cur_node).unwrap();
                let maglie =
                    match self
                        .get_surface()
                        .geodesic_dist(w.coords, w.height_target(), None)
                    {
                        Ok(e) => e,
                        _ => break,
                    } / crochet_consts::get_height();

                outp.push_str(&format!("{}\n", maglie));

                cur_node = nuovo;
            }
            let mut output = File::create(filename).unwrap();
            write!(output, "{}", outp).unwrap();
        }

        fn dot_graph(&self, filename: Option<String>) -> String {
            let outp = format!("{:?}", Dot::new(&self.get_graph()));
            if filename.is_some() {
                let mut output = File::create(filename.unwrap()).unwrap();
                write!(output, "{}", outp).unwrap();
            }
            return outp;
        }

        fn to_string_vertices(&self, filename: Option<String>) -> String {
            let mut my_string = "".to_string();
            let mut cur_node = Some(self.get_start_node());
            loop {
                if cur_node.is_none() {
                    break;
                }
                let p = self
                    .get_surface()
                    .to_3d_point(
                        self.get_graph()
                            .node_weight(cur_node.unwrap())
                            .unwrap()
                            .coords,
                    )
                    .unwrap();
                my_string.push_str(&format!("{} {} {}\n", &p.x, &p.y, &p.z));
                cur_node = self.find_next_course(cur_node.unwrap(), true);
            }

            if filename.is_some() {
                let mut output = File::create(filename.unwrap()).unwrap();
                write!(output, "{}", my_string).unwrap();
            }
            return my_string;
        }

        fn verts_count(&self) -> usize {
            return self.get_graph().node_indices().count();
        }

        fn find_next_course(
            &self,
            start_vert: NodeIndex,
            allow_override: bool,
        ) -> Option<NodeIndex> {
            if allow_override {
                for edge in self.get_graph().edges(start_vert) {
                    if edge.weight().link_type == EdgeType::CourseNextOverrideStart {
                        return Some(edge.target());
                    }
                }
                for edge in self.get_graph().edges(start_vert) {
                    if edge.weight().link_type == EdgeType::CourseNextOverrideEnd {
                        return Some(edge.target());
                    }
                }
            }

            for edge in self.get_graph().edges(start_vert) {
                if edge.weight().link_type == EdgeType::Course {
                    return Some(edge.target());
                }
            }
            return None;
        }

        fn find_next_course_edge(
            &self,
            n: NodeIndex,
            find_start: bool,
            find_end: bool,
        ) -> Option<EdgeIndex> {
            if find_start || find_end {
                for edge in self.get_graph().edges(n) {
                    if (find_start && edge.weight().link_type == EdgeType::CourseNextOverrideStart)
                        || (find_end && edge.weight().link_type == EdgeType::CourseNextOverrideEnd)
                    {
                        return Some(self.get_graph().find_edge(n, edge.target()).unwrap());
                    }
                }
            }

            for edge in self.get_graph().edges(n) {
                if edge.weight().link_type == EdgeType::Course {
                    return Some(self.get_graph().find_edge(n, edge.target()).unwrap());
                }
            }
            return None;
        }

        fn find_prev_course_edge(
            &self,
            n: NodeIndex,
            find_start: bool,
            find_end: bool,
        ) -> Option<EdgeIndex> {
            let edges_to = self
                .get_graph()
                .edges_directed(n, petgraph::Direction::Incoming);
            for e in edges_to {
                if (find_start && e.weight().link_type == EdgeType::CourseNextOverrideStart)
                    || (find_end && e.weight().link_type == EdgeType::CourseNextOverrideEnd)
                {
                    return Some(self.get_graph().find_edge(e.source(), n).unwrap());
                }
            }

            let edges_to = self
                .get_graph()
                .edges_directed(n, petgraph::Direction::Incoming);
            for e in edges_to {
                if e.weight().link_type == EdgeType::Course {
                    return Some(self.get_graph().find_edge(e.source(), n).unwrap());
                }
            }
            return None;
        }

        fn filter_edges(&self, start_vert: NodeIndex, edge_type: EdgeType) -> Vec<NodeIndex> {
            let mut res: Vec<NodeIndex> = Vec::new();
            for edge in self.get_graph().edges(start_vert) {
                if edge.weight().link_type == edge_type {
                    res.push(edge.target());
                }
            }
            return res;
        }

        fn find_all_wales(&self, start_vert: NodeIndex) -> Vec<NodeIndex> {
            let mut dest: Vec<NodeIndex> = Vec::new();
            for edge in self.get_graph().edges(start_vert) {
                if edge.weight().link_type == EdgeType::Wale {
                    dest.push(edge.target());
                }
            }

            return dest;
        }

        fn find_next_wale(&self, start_vert: NodeIndex) -> Option<NodeIndex> {
            let mut wale_edges: Vec<petgraph::graph::EdgeReference<CrochetEdge>> = Vec::new();
            for edge in self.get_graph().edges(start_vert) {
                if edge.weight().link_type == EdgeType::Wale {
                    wale_edges.push(edge);
                }
            }
            // Ordino
            if wale_edges.len() == 0 {
                return None;
            } else if wale_edges.len() == 1 {
                return Some(wale_edges[0].target());
            } else {
                // Quale viene prima o dopo?
                if self
                    .get_graph()
                    .contains_edge(wale_edges[0].target(), wale_edges[1].target())
                {
                    return Some(wale_edges[1].target());
                } else {
                    return Some(wale_edges[0].target());
                }
            }
        }

        /// Restituisce o il wale sovrastante, se è uno, oppure una coppia
        /// ordinata
        ///
        /// Esplode se sono 3.
        fn ordered_wales(&self, start_vert: NodeIndex) -> Vec<NodeIndex> {
            let w = self.find_all_wales(start_vert);
            if w.len() > 2 {
                println!(
                    "Errore a {}",
                    self.get_surface()
                        .to_3d_point(self.node_taucoords(start_vert))
                        .unwrap()
                );
                panic!("Più di 2 wales")
            } else if w.len() <= 1 {
                return w;
            }

            if self.get_graph().find_edge(w[0], w[1]).is_some() {
                return w;
            } else {
                return vec![w[1], w[0]];
            }
        }

        /// Per ogni lato, qual è l'errore percentuale?
        ///
        /// Restituisce un vettore di errori percentuali.
        /// Utile per visualizzazione come istogramma
        fn edge_error_distribution(
            &self,
            filename: Option<String>,
            tipo: Option<EdgeType>,
        ) -> Vec<f64> {
            let mut out_data: Vec<f64> = Vec::new();
            let mut my_string = "".to_string();
            for edge in self.get_graph().edge_references() {
                if let Some(ref tipo_estratto) = tipo {
                    if edge.weight().link_type != tipo_estratto.clone() {
                        continue;
                    }
                }
                let coord_0: TauCoords =
                    self.get_graph().node_weight(edge.source()).unwrap().coords;
                let coord_1: TauCoords =
                    self.get_graph().node_weight(edge.target()).unwrap().coords;

                let target_length: f64 = match edge.weight().link_type {
                    EdgeType::Course => crochet_consts::get_width(),
                    EdgeType::Wale => crochet_consts::get_height(),
                    _ => continue,
                };
                let actual_length: f64 = self
                    .get_surface()
                    .geodesic_dist(coord_0, coord_1, None)
                    .unwrap();

                let percentual_error: f64 = (actual_length - target_length) / target_length;
                out_data.push(percentual_error);

                my_string.push_str(&format!("{}\n", &percentual_error).as_str());
            }

            if filename.is_some() {
                let mut output = File::create(filename.unwrap()).unwrap();
                write!(output, "{}", my_string).unwrap();
            }
            return out_data;
        }

        /// Shear, traslazione laterale rispetto a base
        ///
        /// Quanto si spostala lateramente rispetto
        fn metrics_skew(&self, filename: Option<String>) -> Vec<f64> {
            let mut my_string = "".to_string();
            let mut ret: Vec<f64> = Vec::new();

            for e in self.get_graph().edge_references() {
                if vec![
                    EdgeType::Course,
                    EdgeType::CourseNextOverrideStart,
                    EdgeType::CourseNextOverrideEnd,
                ]
                .contains(&e.weight().link_type)
                {
                    let basi = [e.source(), e.target()];

                    let upnodes = [self.ordered_wales(basi[0]), self.ordered_wales(basi[1])];

                    if upnodes[0].len() == 0 || upnodes[1].len() == 0 {
                        continue;
                    }
                    let v = internal_compute_skew(
                        [
                            self.get_surface()
                                .to_3d_point(self.node_taucoords(basi[0]))
                                .unwrap(),
                            self.get_surface()
                                .to_3d_point(self.node_taucoords(basi[1]))
                                .unwrap(),
                        ],
                        [
                            self.get_surface()
                                .to_3d_point(
                                    self.node_taucoords(upnodes[0].last().unwrap().clone())
                                )
                                .unwrap(),
                            self.get_surface()
                                .to_3d_point(
                                    self.node_taucoords(upnodes[1].first().unwrap().clone())
                                )
                                .unwrap(),
                        ]
                    );

                    ret.push(v);
                    my_string.push_str(&format!(
                        "{}\n",
                        v
                    ));
                }
            }

            if filename.is_some() {
                let mut output = File::create(filename.unwrap()).unwrap();
                write!(output, "{}", my_string).unwrap();
            }


            return ret;
        }

        fn metrics_knitting_direction(&self, filename: Option<String>) {
            fn compute_direction_amount(
                surf: &Box<dyn UVSurface>,
                base_node: TauCoords,
                upnode: TauCoords,
            ) -> f64 {
                let global_knit_direction = surf
                    .get_curve()
                    .first_derivative_normalized(surf.tau_to_u(base_node.tau).unwrap())
                    .unwrap();
                let wale_direction =
                    surf.to_3d_point(upnode).unwrap() - surf.to_3d_point(base_node).unwrap();

                let normal = surf.normal_versor(base_node);

                let knit_dir_projected = (global_knit_direction
                    - global_knit_direction.dot(&normal) * normal)
                    .normalize();
                let wale_projected =
                    (wale_direction - wale_direction.dot(&normal) * normal).normalize();

                return 1.0 - wale_projected.dot(&knit_dir_projected);
            }

            let mut my_string = "".to_string();

            for e in self.get_graph().edge_references() {
                if e.weight().link_type == EdgeType::Wale {
                    my_string.push_str(&format!(
                        "{}\n",
                        compute_direction_amount(
                            self.get_surface(),
                            self.node_taucoords(e.source()),
                            self.node_taucoords(e.target())
                        )
                    ));
                }
            }

            if filename.is_some() {
                let mut output = File::create(filename.unwrap()).unwrap();
                write!(output, "{}", my_string).unwrap();
            }
        }

        /// Parto da start_node e salgo.
        /// Quali nodi segnano l'inizio di ogni riga?
        fn start_rows_nodes(&self, filename: Option<String>) -> String {
            let mut nodes: Vec<NodeIndex> = vec![self.get_start_node()];
            while let Some(next) = self.find_next_wale(nodes.last().unwrap().clone()) {
                nodes.push(next);
            }
            let mut my_string = "".to_string();

            for cur_node in nodes {
                let p = self
                    .get_surface()
                    .to_3d_point(self.get_graph().node_weight(cur_node).unwrap().coords)
                    .unwrap();
                my_string.push_str(&format!("{} {} {}\n", &p.x, &p.y, &p.z));
            }

            if filename.is_some() {
                let mut output = File::create(filename.unwrap()).unwrap();
                write!(output, "{}", my_string).unwrap();
            }
            return my_string;
        }

        // -----------------------------------
        //          DATA TEST SECTION
        // -----------------------------------
        fn print_node(&self, n: NodeIndex) -> String {
            format!("{:?}", self.get_graph().node_weight(n).unwrap())
        }

        fn debug_shortrow_level_nodes(&self, filename: Option<String>, level: i32) -> String {
            let mut my_string = "".to_string();
            for n in self.get_graph().node_indices() {
                if self.get_graph().node_weight(n).unwrap().short_row_level == level {
                    let p = self
                        .get_surface()
                        .to_3d_point(self.get_graph().node_weight(n).unwrap().coords)
                        .unwrap();
                    my_string.push_str(&format!("{} {} {}\n", &p.x, &p.y, &p.z));
                }
            }

            if filename.is_some() {
                let mut output = File::create(filename.unwrap()).unwrap();
                write!(output, "{}", my_string).unwrap();
            }
            return my_string;
        }

        fn edge_middle_point(&self, e: EdgeIndex) -> Point3<f64> {
            let edge_data = self.get_graph().edge_endpoints(e).unwrap();
            let coord_0 = self.get_graph().node_weight(edge_data.0).unwrap();
            let coord_1 = self.get_graph().node_weight(edge_data.1).unwrap();

            return self
                .get_surface()
                .to_3d_point(coord_0.coords.lerp(&coord_1.coords, 0.5))
                .unwrap();
        }

        fn node_taucoords(&self, n: NodeIndex) -> TauCoords {
            self.get_graph().node_weight(n).unwrap().coords
        }

        /// TEST validità grafo
        ///
        /// Da lanciare subito dopo che è stato generato il grafo.
        ///
        /// Opera vari controlli. I controlli sono lanciati solo se ci si trova
        /// in una build di debug per non appesantire il runtime.
        fn test_graph_validity(&self) {
            #[cfg(debug_assertions)]
            {
                // Non ci devono essere archi che puntano a sé stessi
                for ed in self.get_graph().edge_references() {
                    assert_ne!(ed.source(), ed.target());
                }

                // Non ci devono essere nodi con più di due wale in uscita
                // for l in self.get_graph().node_indices() {
                //     assert!(
                //         self.find_all_wales(l).len() <= 2,
                //         "Più di due wale in uscita da {:?}",
                //         self.get_graph().node_weight(l).unwrap().coords
                //     );
                // }

                // Il numero si start_override e end_override dovrebbe coincidere
                let mut count_start: u32 = 0;
                let mut count_end: u32 = 0;
                for e in self.get_graph().edge_weights() {
                    match e.link_type {
                        EdgeType::CourseNextOverrideStart => count_start += 1,
                        EdgeType::CourseNextOverrideEnd => count_end += 1,
                        _ => {}
                    };
                }
                assert_eq!(count_start, count_end);

                // Non devo avere nodi con wale scollegati a parte il primo giro
                for n in self.get_graph().node_indices() {
                    // if self.get_graph().node_weight(n).unwrap().coords.tau < 1.0
                    // ||
                    // self.find_all_wales(n).len() > 0
                    // {
                    //     println!("Errore a {:?}", self.print_node(n));
                    // }
                    let tau: f64 = self.get_graph().node_weight(n).unwrap().coords.tau;
                    let mut found_edges: u32 = 0;
                    for e in self
                        .get_graph()
                        .edges_directed(n, petgraph::Direction::Incoming)
                    {
                        if e.weight().link_type == EdgeType::Wale {
                            found_edges += 1;
                        }
                    }

                    assert!(
                        tau < 1.0 || found_edges > 0,
                        "Errore a {:?}",
                        self.print_node(n)
                    );
                }
            }
        }
    }

    fn static_find_next_course(
        gr: &Graph<CrochetVertex, CrochetEdge>,
        startnode: NodeIndex,
    ) -> Option<NodeIndex> {
        for edge in gr.edges(startnode) {
            if edge.weight().link_type == EdgeType::CourseNextOverrideStart {
                return Some(edge.target());
            } else if edge.weight().link_type == EdgeType::CourseNextOverrideEnd {
                return Some(edge.target());
            }
        }
        for edge in gr.edges(startnode) {
            if edge.weight().link_type == EdgeType::Course {
                return Some(edge.target());
            }
        }
        return None;
    }

    pub struct CapunamanPaperGraphCompliant {
        graph: Graph<CrochetVertex, CrochetEdge>,
        surface: Box<dyn UVSurface>,
        start_node: NodeIndex,
    }
    impl CapunamanPaperGraphCompliant {
        pub fn new(surf: Box<dyn UVSurface>) -> Self {
            let mut gr: Graph<CrochetVertex, CrochetEdge> = Graph::new();
            let n = gr.add_node(CrochetVertex {
                coords: TauCoords { tau: 0.0, v: 0.0 },
                short_row_level: 0,
                target_tau: Some(1.0),
            });

            CapunamanPaperGraphCompliant {
                graph: gr,
                surface: surf,
                start_node: n,
            }
        }
    }
    impl KnittableGraph for CapunamanPaperGraphCompliant {
        fn generate_graph(&mut self) {
            let mut created_nodes: Vec<NodeIndex> = vec![self.get_start_node()];
            let mut lasttau = 0.0;
            loop {
                lasttau = match self.get_surface().width_increment(
                    self.get_surface().helix(lasttau),
                    crochet_consts::get_width(),
                ) {
                    Ok(e) => e.tau,
                    _ => break,
                };
                if lasttau > self.get_surface().tau_max_value() {
                    break;
                }
                #[cfg(debug_assertions)]
                {
                    println!(
                        "{:.3} su {:.3} - {}",
                        lasttau,
                        self.get_surface().tau_max_value(),
                        self.graph.node_count()
                    );
                }

                let h = self.get_surface().helix(lasttau);
                created_nodes.push(self.graph.add_node(CrochetVertex {
                    coords: h,
                    target_tau: Some(h.tau + 1.0),
                    short_row_level: 0,
                }));
            }
            for n in 0..created_nodes.len() - 1 {
                self.graph.add_edge(
                    created_nodes[n],
                    created_nodes[n + 1],
                    CrochetEdge {
                        link_type: EdgeType::Course,
                    },
                );
            }

            let mut added_up: HashSet<NodeIndex> = HashSet::new();
            for nodo_base in created_nodes.iter() {
                let base_coords = self
                    .get_graph()
                    .node_weight(nodo_base.clone())
                    .unwrap()
                    .coords;
                let possibili_wale: Vec<NodeIndex> = self
                    .graph
                    .node_references()
                    .filter(|n| {
                        n.1.coords.tau > base_coords.tau + 0.8
                            && n.1.coords.tau < base_coords.tau + 1.2
                    })
                    .map(|n| n.0)
                    .collect();

                let mut nearest_node: Option<NodeIndex> = None;
                let mut min_distance: Option<f64> = None;
                for next_cycle_node in possibili_wale.iter() {
                    let dist: f64 = self
                        .get_surface()
                        .geodesic_dist(
                            self.get_graph()
                                .node_weight(nodo_base.clone())
                                .unwrap()
                                .coords,
                            self.get_graph()
                                .node_weight(*next_cycle_node)
                                .unwrap()
                                .coords,
                            None,
                        )
                        .unwrap();

                    if min_distance.is_none() || dist < min_distance.unwrap() {
                        min_distance = Some(dist);
                        nearest_node = Some(*next_cycle_node);
                    }
                }

                if nearest_node.is_none() {
                    continue;
                }

                if !self
                    .graph
                    .contains_edge(nodo_base.clone(), nearest_node.unwrap())
                {
                    added_up.insert(nearest_node.unwrap());
                    self.graph.add_edge(
                        nodo_base.clone(),
                        nearest_node.unwrap(),
                        CrochetEdge {
                            link_type: EdgeType::Wale,
                        },
                    );
                }
            }

            // Ci sono nodi in alto a cui non ho collegato niente
            // Scansiono di nuovo quelli base per arrivarci
            // CONTRARIO
            let not_connected = created_nodes.iter().filter(|e| !added_up.contains(e));
            for next_cycle_node in not_connected {
                let mut nearest_node: Option<NodeIndex> = None;
                let mut min_distance: Option<f64> = None;

                let next_c_node_coord = self
                    .get_graph()
                    .node_weight(next_cycle_node.clone())
                    .unwrap()
                    .coords;

                for past_cycle_node in created_nodes
                    .iter()
                    .filter(|n_id| {
                        let c = self
                            .get_graph()
                            .node_weight(*n_id.clone())
                            .unwrap()
                            .coords
                            .tau;
                        c < next_c_node_coord.tau && c > next_c_node_coord.tau - 2.0
                    })
                    .filter(|n| *n != next_cycle_node)
                    .filter(|n| self.find_all_wales(*n.clone()).len() < 2)
                {
                    let dist: f64 = self
                        .get_surface()
                        .geodesic_dist(
                            self.get_graph()
                                .node_weight(*past_cycle_node)
                                .unwrap()
                                .coords,
                            self.get_graph()
                                .node_weight(*next_cycle_node)
                                .unwrap()
                                .coords,
                            None,
                        )
                        .unwrap();

                    if min_distance.is_none() || dist < min_distance.unwrap() {
                        min_distance = Some(dist);
                        nearest_node = Some(*past_cycle_node);
                    }
                }

                if let Some(found) = nearest_node {
                    if !self.graph.contains_edge(found.clone(), *next_cycle_node) {
                        self.graph.add_edge(
                            found,
                            *next_cycle_node,
                            CrochetEdge {
                                link_type: EdgeType::Wale,
                            },
                        );
                    }
                }
            }
            self.test_graph_validity();
        }
        fn get_start_node(&self) -> NodeIndex {
            return self.start_node.clone();
        }
        fn get_graph<'a>(&'a self) -> &'a Graph<CrochetVertex, CrochetEdge> {
            &self.graph
        }
        fn get_surface<'a>(&'a self) -> &'a Box<dyn UVSurface> {
            &self.surface
        }
    }

    pub struct CapunamanPaperGraph {
        graph: Graph<CrochetVertex, CrochetEdge>,
        surface: Box<dyn UVSurface>,
        start_node: NodeIndex,
    }
    impl CapunamanPaperGraph {
        pub fn new(surf: Box<dyn UVSurface>) -> Self {
            let mut gr: Graph<CrochetVertex, CrochetEdge> = Graph::new();
            let n = gr.add_node(CrochetVertex {
                coords: TauCoords { tau: 0.0, v: 0.0 },
                short_row_level: 0,
                target_tau: Some(1.0),
            });

            CapunamanPaperGraph {
                graph: gr,
                surface: surf,
                start_node: n,
            }
        }
        fn fix_graph(&mut self) {
            // FIX per diminuzione subito seguita da aumento.
            // Le trasformo in maglie semplici
            {
                for nodo in self.get_graph().node_indices() {
                    if self
                        .get_graph()
                        .edges_directed(nodo, Direction::Outgoing)
                        .filter(|e| e.weight().link_type == EdgeType::Wale)
                        .count()
                        > 1
                    {
                        // Se mi trovo in un aumento

                        if let Some(edge_2) = self.find_prev_course_edge(nodo, false, false) {
                            if let Some(edge_1) = self.find_prev_course_edge(
                                self.get_graph().edge_endpoints(edge_2).unwrap().0,
                                false,
                                false,
                            ) {
                                let punti: [NodeIndex; 3] = [
                                    self.get_graph().edge_endpoints(edge_1).unwrap().0,
                                    self.get_graph().edge_endpoints(edge_2).unwrap().0,
                                    self.get_graph().edge_endpoints(edge_2).unwrap().1,
                                ];

                                if self.find_next_wale(punti[0]) == self.find_next_wale(punti[1])
                                    && self.find_next_wale(punti[0]).is_some()
                                {
                                    // quella prima è una diminuzione
                                    // tolgo edge non necessari e ne aggiungo uno utile

                                    self.graph.remove_edge(
                                        self.get_graph()
                                            .find_edge(
                                                punti[1],
                                                self.find_next_wale(punti[1]).unwrap(),
                                            )
                                            .unwrap(),
                                    );
                                    let upnode_tmp =
                                        self.ordered_wales(punti[2]).first().unwrap().clone();
                                    self.graph.remove_edge(
                                        self.get_graph().find_edge(punti[2], upnode_tmp).unwrap(),
                                    );

                                    self.graph.add_edge(
                                        punti[1],
                                        upnode_tmp,
                                        CrochetEdge {
                                            link_type: EdgeType::Wale,
                                        },
                                    );
                                }
                            }
                        }
                    }
                }
            }

            // fix
            //   \/\/ → / /
            {
                for nodo in self.get_graph().node_indices() {
                    let next_course: NodeIndex = match self.find_next_course(nodo, false) {
                        Some(e) => e,
                        None => continue,
                    };

                    let base_nodes = [nodo, next_course];
                    let up_prev = self.ordered_wales(base_nodes[0]);
                    let up_next = self.ordered_wales(base_nodes[1]);

                    if up_prev.len() != 2 || up_next.len() != 2 {
                        continue;
                    }

                    if up_prev[1] == up_next[0] {
                        self.graph.remove_edge(
                            self.get_graph()
                                .find_edge(base_nodes[1], up_next[0])
                                .unwrap(),
                        );
                    }
                }
            }

            // doppioni wale?
            {
                for nodo in self.get_graph().node_indices() {
                    let original_wales = self.find_all_wales(nodo);
                    let mut wales = original_wales.to_vec();

                    wales.sort();

                    for i in 1..wales.len() {
                        if wales[i] == wales[i - 1] {
                            self.graph
                                .remove_edge(self.get_graph().find_edge(nodo, wales[i]).unwrap());
                            break;
                        }
                    }
                }
            }

            // fix
            //   |X → | |
            // {
            //     for nodo in self.get_graph().node_indices()
            //     {
            //         if self.get_graph().edges_directed(
            //             nodo, Direction::Outgoing
            //         ).filter(|e| e.weight().link_type == EdgeType::Wale).count()>1
            //         {
            //             // Mi trovo nel primo aumento

            //         }
            //     }
            // }
        }
    }
    impl KnittableGraph for CapunamanPaperGraph {
        fn generate_graph(&mut self) {
            // let mut last_added_node: NodeIndex = self.get_start_node();
            let mut created_nodes: Vec<Vec<NodeIndex>> = vec![vec![self.get_start_node()]];
            for i in 1..self.get_surface().tau_max_value().floor() as usize
            {
                created_nodes.push(Vec::new());
            }

            let mut last_generated_cycle: Vec<NodeIndex> = vec![];
            for i in 0..created_nodes.len() {
                #[cfg(debug_assertions)]
                {
                    println!("{}", i);
                }
                let mut current_generated_cycle: Vec<NodeIndex> = Vec::new();
                // Genero nuovi nodi nel ciclo attuale
                // Prima di tutto calcolo la distanza
                let cycle_length: f64 = self
                    .get_surface()
                    .geodesic_dist(
                        TauCoords {
                            tau: i as f64,
                            v: 0.0,
                        },
                        TauCoords {
                            tau: i as f64 + 0.25,
                            v: 0.25,
                        },
                        Some(20),
                    )
                    .unwrap()
                    + self
                        .get_surface()
                        .geodesic_dist(
                            TauCoords {
                                tau: i as f64 + 0.25,
                                v: 0.25,
                            },
                            TauCoords {
                                tau: i as f64 + 0.5,
                                v: 0.5,
                            },
                            Some(20),
                        )
                        .unwrap()
                    + self
                        .get_surface()
                        .geodesic_dist(
                            TauCoords {
                                tau: i as f64 + 0.5,
                                v: 0.5,
                            },
                            TauCoords {
                                tau: i as f64 + 0.75,
                                v: 0.75,
                            },
                            Some(20),
                        )
                        .unwrap()
                    + self
                        .get_surface()
                        .geodesic_dist(
                            TauCoords {
                                tau: i as f64 + 0.75,
                                v: 0.75,
                            },
                            TauCoords {
                                tau: i as f64 + 1.0,
                                v: 0.0,
                            },
                            Some(20),
                        )
                        .unwrap();

                let n_stitches: i32 = (cycle_length / crochet_consts::get_width()).round() as i32;

                for j in 0..n_stitches
                {
                    let cur_v = j as f64 / n_stitches as f64;
                    if i!=0 || j!=0
                    {
                        created_nodes[i].push(
                            self.graph.add_node(CrochetVertex
                            {
                                coords: TauCoords{tau: i as f64 + cur_v, v: cur_v},
                                short_row_level: 0,
                                target_tau: Some(i as f64 + 1.0)
                            })
                        );
                    }
                }
            }

            for riga in created_nodes.iter()
            {
                for i in 1..riga.len()
                {
                    self.graph.add_edge(riga[i-1], riga[i],CrochetEdge{link_type: EdgeType::Course});
                }
            }
            for i in 1..created_nodes.len()
            {
                self.graph.add_edge(
                    created_nodes[i-1].last().unwrap().clone(),
                    created_nodes[i].first().unwrap().clone(),
                    CrochetEdge{link_type: EdgeType::Course}
                );
            }


            // Li collego con Wale
            for i in 0..created_nodes.len()-1
            {
                let mut added_up: HashSet<NodeIndex> = HashSet::new();
                for (down_index, nodo_base) in created_nodes[i].iter().enumerate() {
                    let base_coords = self
                        .get_graph()
                        .node_weight(nodo_base.clone())
                        .unwrap()
                        .coords;
                    let possibili_wale: Vec<NodeIndex> = created_nodes[i+1].to_vec();
    
                    let mut nearest_node: Option<NodeIndex> = None;
                    let mut min_distance: Option<f64> = None;
                    for (up_index, next_cycle_node) in possibili_wale.iter().enumerate() {
                        let dist: f64 = self
                            .get_surface()
                            .geodesic_dist(
                                self.get_graph()
                                    .node_weight(nodo_base.clone())
                                    .unwrap()
                                    .coords,
                                self.get_graph()
                                    .node_weight(*next_cycle_node)
                                    .unwrap()
                                    .coords,
                                None,
                            )
                            .unwrap();
                        
                        if down_index as f64 / created_nodes[i].len() as f64 > 0.6
                        &&
                        up_index as f64  / (possibili_wale.len() as f64) < 0.3
                        {
                            continue;
                        }
                        if min_distance.is_none()
                            ||
                            (
                                dist < min_distance.unwrap()
                            )
                        {
                            min_distance = Some(dist);
                            nearest_node = Some(*next_cycle_node);
                        }
                    }
    
                    if nearest_node.is_none() {
                        continue;
                    }
    
                    if !self
                        .graph
                        .contains_edge(nodo_base.clone(), nearest_node.unwrap())
                    {
                        added_up.insert(nearest_node.unwrap());
                        self.graph.add_edge(
                            nodo_base.clone(),
                            nearest_node.unwrap(),
                            CrochetEdge {
                                link_type: EdgeType::Wale,
                            },
                        );
                    }
                }
    
                // Ci sono nodi in alto a cui non ho collegato niente
                // Scansiono di nuovo quelli base per arrivarci
                // CONTRARIO
                let not_connected = created_nodes[i+1].iter().filter(|e| !added_up.contains(e));
                for next_cycle_node in not_connected {
                    let mut nearest_node: Option<NodeIndex> = None;
                    let mut min_distance: Option<f64> = None;
    
                    let next_c_node_coord = self
                        .get_graph()
                        .node_weight(next_cycle_node.clone())
                        .unwrap()
                        .coords;
    
                    for past_cycle_node in created_nodes[i]
                        .iter()
                        .filter(|n| self.find_all_wales(*n.clone()).len() < 2)
                    {
                        let dist: f64 = self
                            .get_surface()
                            .geodesic_dist(
                                self.get_graph()
                                    .node_weight(*past_cycle_node)
                                    .unwrap()
                                    .coords,
                                self.get_graph()
                                    .node_weight(*next_cycle_node)
                                    .unwrap()
                                    .coords,
                                None,
                            )
                            .unwrap();
    
                        if min_distance.is_none() || (
                            dist < min_distance.unwrap()
                            &&
                            self.get_graph()
                            .node_weight(past_cycle_node.clone()).unwrap().coords.v
                            -
                            self.get_graph()
                            .node_weight(*next_cycle_node).unwrap().coords.v
                            < 0.75 
                            ) {
                            min_distance = Some(dist);
                            nearest_node = Some(*past_cycle_node);
                        }
                    }
    
                    if let Some(found) = nearest_node {
                        if !self.graph.contains_edge(found.clone(), *next_cycle_node)
                        &&
                        self.find_all_wales(found.clone()).len()<2
                        {
                            self.graph.add_edge(
                                found,
                                *next_cycle_node,
                                CrochetEdge {
                                    link_type: EdgeType::Wale,
                                },
                            );
                        }
                    }
                }

                // Rimasti nodi senza wale?
                for nodo_base in created_nodes[i].iter()
                {
                    if self.find_all_wales(nodo_base.clone()).len()==0
                    {
                        let possibili_wale: Vec<NodeIndex> = created_nodes[i+1].to_vec();
    
                        let mut nearest_node: Option<NodeIndex> = None;
                        let mut min_distance: Option<f64> = None;
                        for next_cycle_node in possibili_wale.iter() {
                            let dist: f64 = self
                                .get_surface()
                                .geodesic_dist(
                                    self.get_graph()
                                        .node_weight(nodo_base.clone())
                                        .unwrap()
                                        .coords,
                                    self.get_graph()
                                        .node_weight(*next_cycle_node)
                                        .unwrap()
                                        .coords,
                                    None,
                                )
                                .unwrap();
        
                            if min_distance.is_none() || dist < min_distance.unwrap() {
                                min_distance = Some(dist);
                                nearest_node = Some(*next_cycle_node);
                            }
                        }

                        if let Some(ne) = nearest_node
                        {
                            self.graph.add_edge(nodo_base.clone(), ne, CrochetEdge {
                                link_type: EdgeType::Wale,
                            });
                        }
                    }
                }
            }

            self.test_graph_validity();
            self.fix_graph();
        }
        fn get_start_node(&self) -> NodeIndex {
            return self.start_node.clone();
        }
        fn get_graph<'a>(&'a self) -> &'a Graph<CrochetVertex, CrochetEdge> {
            &self.graph
        }
        fn get_surface<'a>(&'a self) -> &'a Box<dyn UVSurface> {
            &self.surface
        }
    }

    pub struct CrochetGraph2 {
        graph: Graph<CrochetVertex, CrochetEdge>,
        surface: Box<dyn UVSurface>,
        start_node: NodeIndex,
    }
    impl CrochetGraph2 {
        pub fn new(surf: Box<dyn UVSurface>) -> Self {
            let mut gr: Graph<CrochetVertex, CrochetEdge> = Graph::new();

            // Primo nodo
            let start_node = gr.add_node(CrochetVertex {
                coords: TauCoords { tau: 0.0, v: 0.0 },
                short_row_level: 0,
                target_tau: Some(1.0),
            });

            // Primo giro
            let cycle_length: f64 = surf
                .geodesic_dist(
                    TauCoords { tau: 0.0, v: 0.0 },
                    TauCoords { tau: 0.25, v: 0.25 },
                    Some(5),
                )
                .unwrap()
                + surf
                    .geodesic_dist(
                        TauCoords { tau: 0.25, v: 0.25 },
                        TauCoords { tau: 0.5, v: 0.5 },
                        Some(5),
                    )
                    .unwrap()
                + surf
                    .geodesic_dist(
                        TauCoords { tau: 0.5, v: 0.5 },
                        TauCoords { tau: 0.75, v: 0.75 },
                        Some(5),
                    )
                    .unwrap()
                + surf
                    .geodesic_dist(
                        TauCoords { tau: 0.75, v: 0.75 },
                        TauCoords { tau: 0.999, v: 0.0 },
                        Some(5),
                    )
                    .unwrap();

            let n_stitches: i32 = (cycle_length / crochet_consts::get_width()).round() as i32;
            assert!(n_stitches > 3);

            let step: f64 = 1.0 / n_stitches as f64;

            let mut prev_added = start_node;

            for j in 1..n_stitches {
                let cur_v = step * j as f64;
                // Il primissimo nodo era già stato aggiunto
                // evitare doppioni in (0, 0)
                let new_node = gr.add_node(CrochetVertex {
                    coords: TauCoords {
                        tau: cur_v,
                        v: cur_v,
                    },
                    short_row_level: 0,
                    target_tau: Some(cur_v + 1.0),
                });
                gr.add_edge(
                    prev_added,
                    new_node,
                    CrochetEdge {
                        link_type: EdgeType::Course,
                    },
                );
                prev_added = new_node;
            }

            CrochetGraph2 {
                graph: gr,
                surface: surf,
                start_node: start_node,
            }
        }

        fn snap(
            surf: &Box<dyn UVSurface>,
            coord: TauCoords,
            target: TauCoords,
        ) -> Result<TauCoords, CrochetError> {
            //const PERCENTAGE_ERROR_ALLOWED: f64 = 0.2;
            let dist = surf.geodesic_dist(coord, target, None)?;
            let limit = crochet_consts::get_height() * crochet_consts::get_snap_treshold();
            if
            // coord.tau > target.tau
            // &&
            dist <= limit {
                return Ok(target);
            } else {
                return Ok(coord);
            }
        }

        fn percentual_threshold_offset(
            surf: &Box<dyn UVSurface>,
            below: [TauCoords; 2],
            current: [TauCoords; 2],
        ) -> Result<f64, CrochetError> {
            // const K: f64 = 2.95;
            // let a = surf.geodesic_dist(below[0], below[1], None)?;
            // let b = surf.geodesic_dist(current[0], current[1], None)?;
            let v1 = surf.to_3d_point(current[0])? - surf.to_3d_point(below[0])?;
            let v2 = surf.to_3d_point(current[1])? - surf.to_3d_point(below[1])?;

            let cross = v1.normalize().cross(&v2.normalize());
            return Ok(-cross.dot(&surf.normal_versor(below[0])));
        }

        fn fix_graph(&mut self) {
            // FIX per diminuzione subito seguita da aumento.
            // Le trasformo in maglie semplici
            {
                for nodo in self.get_graph().node_indices() {
                    if self
                        .get_graph()
                        .edges_directed(nodo, Direction::Outgoing)
                        .filter(|e| e.weight().link_type == EdgeType::Wale)
                        .count()
                        > 1
                    {
                        // Se mi trovo in un aumento

                        if let Some(edge_2) = self.find_prev_course_edge(nodo, false, false) {
                            if let Some(edge_1) = self.find_prev_course_edge(
                                self.get_graph().edge_endpoints(edge_2).unwrap().0,
                                false,
                                false,
                            ) {
                                let punti: [NodeIndex; 3] = [
                                    self.get_graph().edge_endpoints(edge_1).unwrap().0,
                                    self.get_graph().edge_endpoints(edge_2).unwrap().0,
                                    self.get_graph().edge_endpoints(edge_2).unwrap().1,
                                ];

                                if self.find_next_wale(punti[0]) == self.find_next_wale(punti[1])
                                    && self.find_next_wale(punti[0]).is_some()
                                {
                                    // quella prima è una diminuzione
                                    // tolgo edge non necessari e ne aggiungo uno utile

                                    self.graph.remove_edge(
                                        self.get_graph()
                                            .find_edge(
                                                punti[1],
                                                self.find_next_wale(punti[1]).unwrap(),
                                            )
                                            .unwrap(),
                                    );
                                    let upnode_tmp =
                                        self.ordered_wales(punti[2]).first().unwrap().clone();
                                    self.graph.remove_edge(
                                        self.get_graph().find_edge(punti[2], upnode_tmp).unwrap(),
                                    );

                                    self.graph.add_edge(
                                        punti[1],
                                        upnode_tmp,
                                        CrochetEdge {
                                            link_type: EdgeType::Wale,
                                        },
                                    );
                                }
                            }
                        }
                    }
                }
            }

            // fix
            //   \/\/ → / /
            {
                for nodo in self.get_graph().node_indices() {
                    let next_course: NodeIndex = match self.find_next_course(nodo, false) {
                        Some(e) => e,
                        None => continue,
                    };

                    let base_nodes = [nodo, next_course];
                    let up_prev = self.ordered_wales(base_nodes[0]);
                    let up_next = self.ordered_wales(base_nodes[1]);

                    if up_prev.len() != 2 || up_next.len() != 2 {
                        continue;
                    }

                    if up_prev[1] == up_next[0] {
                        self.graph.remove_edge(
                            self.get_graph()
                                .find_edge(base_nodes[1], up_next[0])
                                .unwrap(),
                        );
                    }
                }
            }

            // doppioni wale?
            {
                for nodo in self.get_graph().node_indices() {
                    let original_wales = self.find_all_wales(nodo);
                    let mut wales = original_wales.to_vec();

                    wales.sort();

                    for i in 1..wales.len() {
                        if wales[i] == wales[i - 1] {
                            self.graph
                                .remove_edge(self.get_graph().find_edge(nodo, wales[i]).unwrap());
                            break;
                        }
                    }
                }
            }

            // fix
            //   |X → | |
            // {
            //     for nodo in self.get_graph().node_indices()
            //     {
            //         if self.get_graph().edges_directed(
            //             nodo, Direction::Outgoing
            //         ).filter(|e| e.weight().link_type == EdgeType::Wale).count()>1
            //         {
            //             // Mi trovo nel primo aumento

            //         }
            //     }
            // }
        }

        /// Smussa anche i vertici nelle short row
        fn undeform(&mut self) {
            let gr = self.get_graph();
            let nodes_to_scan: Vec<(
                Vec<NodeIndex>,
                NodeIndex,
                Vec<NodeIndex>,
                Vec<NodeIndex>,
                Vec<NodeIndex>,
            )> = self
                .get_graph()
                .node_references()
                .filter(|n| n.1.short_row_level != 0)
                .map(|e| e.0)
                .map(|e| {
                    (
                        gr.edges_directed(e, Direction::Incoming)
                            .filter(|e| e.weight().link_type == EdgeType::Wale)
                            .map(|e| e.source())
                            .collect(),
                        e,
                        gr.edges_directed(e, Direction::Outgoing)
                            .filter(|e| e.weight().link_type == EdgeType::Wale)
                            .map(|e| e.target())
                            .collect(),
                        gr.edges_directed(e, Direction::Incoming)
                            .filter(|e| {
                                e.weight().link_type == EdgeType::Course
                                    || e.weight().link_type == EdgeType::CourseNextOverrideStart
                                    || e.weight().link_type == EdgeType::CourseNextOverrideEnd
                            })
                            .map(|e| e.source())
                            .collect(),
                        gr.edges_directed(e, Direction::Outgoing)
                            .filter(|e| {
                                e.weight().link_type == EdgeType::Course
                                    || e.weight().link_type == EdgeType::CourseNextOverrideStart
                                    || e.weight().link_type == EdgeType::CourseNextOverrideEnd
                            })
                            .map(|e| e.target())
                            .collect(),
                    )
                })
                .collect();

            // Lerp dei due limiti rispetto ai nodi sopra e sotto
            {
                for l in nodes_to_scan.iter().clone() {
                    let nodi_sotto = &l.0;
                    let nodi_sopra = &l.2;
                    let nodi_prev = &l.3;
                    let nodi_next = &l.4;
                    if nodi_sotto.len() == 0 || nodi_sopra.len() == 0 {
                        continue;
                    }

                    let avg_down_coord: TauCoords = match nodi_sotto.len() {
                        1 => self.get_graph().node_weight(nodi_sotto[0]).unwrap().coords,
                        2 => self
                            .get_graph()
                            .node_weight(nodi_sotto[0])
                            .unwrap()
                            .coords
                            .lerp(
                                &self.get_graph().node_weight(nodi_sotto[1]).unwrap().coords,
                                0.5,
                            ),
                        _ => self.get_graph().node_weight(nodi_sotto[0]).unwrap().coords,
                    };
                    let avg_up_coord: TauCoords = match nodi_sopra.len() {
                        1 => self.get_graph().node_weight(nodi_sopra[0]).unwrap().coords,
                        2 => self
                            .get_graph()
                            .node_weight(nodi_sopra[0])
                            .unwrap()
                            .coords
                            .lerp(
                                &self.get_graph().node_weight(nodi_sopra[1]).unwrap().coords,
                                0.5,
                            ),
                        _ => self.get_graph().node_weight(nodi_sopra[0]).unwrap().coords,
                    };
                    let avg_prev_coord: Option<TauCoords> = match nodi_prev.len() {
                        0 => None,
                        1 => Some(self.get_graph().node_weight(nodi_prev[0]).unwrap().coords),
                        2 => Some(
                            self.get_graph()
                                .node_weight(nodi_prev[0])
                                .unwrap()
                                .coords
                                .lerp(
                                    &self.get_graph().node_weight(nodi_prev[1]).unwrap().coords,
                                    0.5,
                                ),
                        ),
                        _ => Some(self.get_graph().node_weight(nodi_prev[0]).unwrap().coords),
                    };
                    let avg_next_coord: Option<TauCoords> = match nodi_next.len() {
                        0 => None,
                        1 => Some(self.get_graph().node_weight(nodi_next[0]).unwrap().coords),
                        2 => Some(
                            self.get_graph()
                                .node_weight(nodi_next[0])
                                .unwrap()
                                .coords
                                .lerp(
                                    &self.get_graph().node_weight(nodi_next[1]).unwrap().coords,
                                    0.5,
                                ),
                        ),
                        _ => Some(self.get_graph().node_weight(nodi_next[0]).unwrap().coords),
                    };

                    let wale_avg = avg_down_coord.lerp(&avg_up_coord, 0.5);
                    if avg_prev_coord.is_none() || avg_next_coord.is_none() {
                        self.graph.node_weight_mut(l.1).unwrap().coords = wale_avg;
                    } else {
                        let horiz_avg = avg_prev_coord.unwrap().lerp(&avg_next_coord.unwrap(), 0.5);
                        self.graph.node_weight_mut(l.1).unwrap().coords =
                            wale_avg.lerp(&horiz_avg, 0.5);
                    }
                }
            }

            // Correzione estremi short row
            for (start, end) in self.shortrows_limits() {
                for n in [start, end] {
                    let sotto: Vec<NodeIndex> = self
                        .get_graph()
                        .edges_directed(n, Direction::Incoming)
                        .filter(|e| e.weight().link_type == EdgeType::Wale)
                        .map(|e| e.source())
                        .collect();

                    let sopra: Vec<NodeIndex> = self
                        .get_graph()
                        .edges_directed(n, Direction::Outgoing)
                        .filter(|e| e.weight().link_type == EdgeType::Wale)
                        .map(|e| e.target())
                        .collect();

                    if sotto.len() == 0 || sopra.len() == 0 {
                        continue;
                    }

                    let avg_down_coord: TauCoords = match sotto.len() {
                        1 => self.get_graph().node_weight(sotto[0]).unwrap().coords,
                        2 => self
                            .get_graph()
                            .node_weight(sotto[0])
                            .unwrap()
                            .coords
                            .lerp(&self.get_graph().node_weight(sotto[1]).unwrap().coords, 0.5),
                        _ => self.get_graph().node_weight(sotto[0]).unwrap().coords,
                    };
                    let avg_up_coord: TauCoords = match sopra.len() {
                        1 => self.get_graph().node_weight(sopra[0]).unwrap().coords,
                        2 => self
                            .get_graph()
                            .node_weight(sopra[0])
                            .unwrap()
                            .coords
                            .lerp(&self.get_graph().node_weight(sopra[1]).unwrap().coords, 0.5),
                        _ => self.get_graph().node_weight(sopra[0]).unwrap().coords,
                    };

                    let wale_avg = avg_down_coord.lerp(&avg_up_coord, 0.5);
                    self.graph.node_weight_mut(n).unwrap().coords = wale_avg;
                }
            }
        }
    }

    impl KnittableGraph for CrochetGraph2 {
        fn generate_graph(&mut self) {
            /// Quanti nodi devo generare a salire?
            ///
            /// Restituisce TauCoords dei nodi su questa verticale,
            /// compreso quello passato come _start_
            fn nodes_to_generate_height(
                surf: &Box<dyn UVSurface>,
                start: TauCoords,
                end: TauCoords,
                need_to_close_shortrows: bool,
            ) -> Result<Vec<TauCoords>, CrochetError> {
                let mut ret_vec: Vec<TauCoords> = Vec::new();
                ret_vec.push(start);
                let mut last_added = start;
                loop {
                    if need_to_close_shortrows && ret_vec.len() == 2 {
                        return Ok(ret_vec);
                    }
                    let new_coord: TauCoords = CrochetGraph2::snap(
                        surf,
                        match surf.height_increment(last_added, crochet_consts::get_height()) {
                            Ok(e) => e,
                            Err(_) => break,
                        },
                        end,
                    )?;
                    if new_coord.tau > end.tau {
                        break;
                    }
                    // Se i nodi generati sono coincidenti?
                    if surf.geodesic_dist(new_coord, ret_vec.last().unwrap().clone(), None)?
                        < crochet_consts::get_height() * 0.001
                    {
                        break;
                    }

                    ret_vec.push(new_coord);
                    last_added = new_coord;
                }

                if ret_vec.len() < 2 {
                    ret_vec.push(surf.height_increment(ret_vec[0], crochet_consts::get_height())?);
                }

                return Ok(ret_vec);
            }

            fn yarn_relaxation(
                gr: &mut Graph<CrochetVertex, CrochetEdge>,
                surf: &Box<dyn UVSurface>,
                nodes_to_smooth: &mut Vec<NodeIndex>,
            ) {
                if nodes_to_smooth.len() != SMOOTH_WINDOW_SIZE {
                    return;
                }
                let mut left_average: TauCoords =
                    gr.node_weight(nodes_to_smooth[0]).unwrap().coords;
                let mut right_average: TauCoords = gr
                    .node_weight(nodes_to_smooth.last().unwrap().clone())
                    .unwrap()
                    .coords;

                let middle_index = nodes_to_smooth.len() / 2;

                for i in 1..nodes_to_smooth.len() / 2 {
                    left_average =
                        left_average.lerp(&gr.node_weight(nodes_to_smooth[i]).unwrap().coords, 0.5);
                }
                for i in (nodes_to_smooth.len() / 2 + 1..nodes_to_smooth.len() - 1).rev() {
                    right_average = right_average
                        .lerp(&gr.node_weight(nodes_to_smooth[i]).unwrap().coords, 0.5);
                }

                let mut new_middle_coords = left_average.lerp(&right_average, 0.5);

                let current_coords: TauCoords = gr
                    .node_weight(nodes_to_smooth[middle_index])
                    .unwrap()
                    .coords;
                let current_target: f64 = gr
                    .node_weight(nodes_to_smooth[middle_index])
                    .unwrap()
                    .target_tau
                    .unwrap();
                let mut new_coords = current_coords.lerp(&new_middle_coords, 0.5);

                let middle_prev = gr
                    .node_weight(nodes_to_smooth[0])
                    .unwrap()
                    .coords
                    .lerp(&gr.node_weight(nodes_to_smooth[2]).unwrap().coords, 0.5);
                gr.node_weight_mut(nodes_to_smooth[1]).unwrap().coords =
                    middle_prev.lerp(&gr.node_weight(nodes_to_smooth[1]).unwrap().coords, 0.8);

                if let Ok(dist_tmp) = surf.geodesic_dist(
                    current_coords,
                    current_coords + TauCoords { tau: 1.0, v: 0.0 },
                    None,
                ) {
                    if dist_tmp < crochet_consts::get_height() * 0.95
                    // &&
                    // dist_tmp > crochet_consts::get_height()*0.95
                    {
                        if let Some(below_edge) = gr
                            .edges_directed(nodes_to_smooth[middle_index], Direction::Incoming)
                            .filter(|e| e.weight().link_type == EdgeType::Wale)
                            .next()
                        {
                            let below_node_w = gr.node_weight(below_edge.source()).unwrap();
                            let dist: f64 = surf
                                .geodesic_dist(below_node_w.coords, new_coords, None)
                                .unwrap();
                            let delta_tau: f64 = new_coords.tau - below_node_w.coords.tau;

                            let wanted_length = crochet_consts::get_height();

                            if dist < crochet_consts::get_height() {
                                new_coords.tau = (below_node_w.coords.tau
                                    + ((delta_tau / dist) * crochet_consts::get_height()))
                                .clamp(0.0, surf.tau_max_value());
                            }
                        }
                    }
                }

                // let sotto: Vec<NodeIndex> = gr.edges_directed(
                //     nodes_to_smooth[middle_index], Direction::Incoming
                // )
                // .filter(|e|e.weight().link_type == EdgeType::Wale)
                // .map(|e|e.source())
                // .collect();
                // let down_coords: TauCoords;
                // {
                //     if sotto.len()==0
                //     {
                //         down_coords = new_coords;
                //     }
                //     else
                //     {
                //         down_coords = match sotto.len()
                //         {
                //             1 => gr.node_weight(sotto[0]).unwrap().coords,
                //             2 => gr.node_weight(sotto[0]).unwrap().coords.lerp(
                //                 &gr.node_weight(sotto[1]).unwrap().coords,
                //                 0.5
                //             ),
                //             _ => gr.node_weight(sotto[0]).unwrap().coords
                //         };
                //     }
                // }
                // new_coords.v = down_coords.lerp(&new_coords, 0.5).v;

                gr.node_weight_mut(nodes_to_smooth[middle_index])
                    .unwrap()
                    .coords = new_coords;
                gr.node_weight_mut(nodes_to_smooth[middle_index])
                    .unwrap()
                    .target_tau = Some(
                    TauCoords {
                        tau: current_target,
                        v: current_target,
                    }
                    .move_right_helix(new_middle_coords.v - current_coords.v)
                    .tau,
                );
            }

            fn should_split(
                surf: &Box<dyn UVSurface>,
                column_start: [TauCoords; 2],
                column_middle: [TauCoords; 2],
                column_end: [TauCoords; 2],
            ) -> bool {
                let coords_start = column_start[1];
                let coords_end = column_end[1];

                let middle_virtual = coords_start.lerp(&coords_end, 0.5);
                let dist = surf
                    .geodesic_dist(coords_start, middle_virtual, None)
                    .unwrap();

                let perc_threshold = CrochetGraph2::percentual_threshold_offset(
                    surf,
                    [column_start[0], column_middle[0]],
                    [column_start[1], column_middle[1]],
                )
                .unwrap();
                // return dist > crochet_consts::target_to_split_width();
                let l_normalized = dist / crochet_consts::get_width();
                return l_normalized
                + crochet_consts::dithering(coords_start.v)
                >
                crochet_consts::target_to_split_width(perc_threshold);
            }

            fn should_join(
                surf: &Box<dyn UVSurface>,
                column_start: [TauCoords; 2],
                column_middle: [TauCoords; 2],
                column_end: [TauCoords; 2]
            ) -> bool {
                let coords_start = column_start[1];
                let coords_middle = column_middle[1];

                let mut new_coord = coords_start;
                new_coord.tau = coords_middle.tau;

                let dist = surf.geodesic_dist(column_start[1], column_middle[1], None).unwrap();

                let perc_threshold = CrochetGraph2::percentual_threshold_offset(
                    surf,
                    [column_start[0], column_middle[0]],
                    [column_start[1], column_middle[1]],
                )
                .unwrap();

                let l_normalized = dist / crochet_consts::get_width();

                let very_close = l_normalized
                + crochet_consts::dithering(coords_start.v)
                < crochet_consts::target_to_join_width(perc_threshold);
                
                if very_close {
                    return true;
                }

                let delta_v = coords_middle.v - coords_start.v;
                if (delta_v >= 0.0 && delta_v < 0.5) || (delta_v < -0.5) {
                    return false;
                }
                return true;
            }

            /// genera le coordinate a cui dovrò creare dei nodi
            fn generate_current_step_coords(
                surf: &Box<dyn UVSurface>,
                column_left: &Vec<TauCoords>,
                column_middle: &Vec<TauCoords>,
                column_right: &Vec<TauCoords>,
                past_was_joined: usize, //Il nodo precedente era stato joinato?
                past_was_splitted: usize,
            ) -> (Vec<Vec<TauCoords>>, usize, usize) //Primo indice altezza, secondo indice orizzontale
            {
                let mut ret: Vec<Vec<TauCoords>> = Vec::new();
                ret.push(Vec::new());

                ret[0].push(column_middle[0]);

                if past_was_joined == 0
                    && past_was_splitted == 0
                    && should_join(
                        surf,
                        [column_left[0], column_left[1]],
                        [column_middle[0], column_middle[1]],
                        [column_right[0], column_right[1]]
                    )
                    && (column_left.len() == column_middle.len()
                        && column_middle.len() == column_right.len())
                {
                    return (ret, 2, 0);
                }


                let mut split_done = false;
                for current_level in 1..column_middle.len() {
                    // Devo aggiungere qualcosa?

                    // aggiungo la singola nuova coordinata
                    ret.push(vec![column_middle[current_level]]);

                    // Quelli ai lati esistono?
                    if column_middle.len() == column_left.len()
                        && column_middle.len() == column_right.len()
                    {
                        //Divido in più maglie?
                        // Coordinate dei nodi che dovrò aggiungere ora
                        let mut coords = vec![
                            column_left[current_level],
                            column_middle[current_level],
                            column_right[current_level],
                        ];
                        loop {
                            let mut splitted_one = false;
                            let mut j: usize = 1;
                            while j < coords.len() - 1 {
                                if past_was_joined == 0
                                &&
                                past_was_splitted == 0
                                    && should_split(
                                        surf,
                                        [column_left[current_level - 1], coords[j - 1]],
                                        [column_middle[current_level - 1], coords[j]],
                                        [column_right[current_level - 1], coords[j + 1]],
                                    )
                                {
                                    splitted_one = true;
                                    split_done = true;
                                    let c1 = coords[j - 1].lerp(&coords[j + 1], 0.33333);
                                    let c2 = coords[j - 1].lerp(&coords[j + 1], 0.66666);

                                    coords[j] = c1;
                                    coords.insert(j + 1, c2);
                                    j += 1;
                                }
                                j += 1;
                            }
                            if !splitted_one {
                                break;
                            }
                        }

                        ret[current_level].pop();
                        for x in 1..coords.len() - 1 {
                            ret[current_level].push(coords[x]);
                        }
                    } else {
                        // Per ora lascio la sola maglia qua
                        // Dovrò farne una per ognuna del livello precedente
                        if ret[current_level - 1].len() > 1 {
                            let delta_tau = column_middle[current_level].tau
                                - column_middle[current_level - 1].tau;
                            ret[current_level] = Vec::new();

                            for j in 0..ret[current_level - 1].len() {
                                let new_value = ret[current_level - 1][j]
                                    + TauCoords {
                                        tau: delta_tau,
                                        v: 0.0,
                                    };
                                ret[current_level].push(new_value);
                            }
                        }
                    }
                }

                for i in 0..ret.len() {
                    while ret[i].len() > 2 {
                        ret[i].pop();
                    }
                }
                let was_splitted = ret[1].len() > 1;

                let p_join: usize = match past_was_joined {
                    2 => 1,
                    1 => 0,
                    _ => 0
                };

                let past_was_splitted = if split_done {2} else {0};
                let p_split: usize = match past_was_splitted {
                    2 => 1,
                    1 => 0,
                    _ => 0
                };

                return (ret, p_join, p_split);
            }

            /// Dalle coordinate, genera i nodi
            /// e li collega tra loro. Forma un "patch"
            ///
            /// Restituisce i nodi al margine destro e sinistro
            ///
            /// Il secondo campo della tupla indica i nodi in cima
            /// che andranno smussati
            fn generate_nodes(
                gr: &mut Graph<CrochetVertex, CrochetEdge>,
                surf: &Box<dyn UVSurface>,
                dati_coords: &Vec<Vec<TauCoords>>,
                base_node: NodeIndex,
                previous_column_nodes: &Vec<NodeIndex>,
            ) -> ([Vec<NodeIndex>; 2], Vec<NodeIndex>) {
                let mut built_nodes: Vec<Vec<NodeIndex>> = vec![vec![base_node]];
                let base_weight: CrochetVertex = gr.node_weight(base_node).unwrap().clone();
                let base_target: TauCoords = gr.node_weight(base_node).unwrap().height_target();

                if dati_coords.len() == 1 {
                    let medium_v: TauCoords = gr
                        .node_weight(previous_column_nodes[0])
                        .unwrap()
                        .coords
                        .lerp(&gr.node_weight(base_node).unwrap().coords, 0.5);
                    let mut nodes: Vec<NodeIndex> = vec![base_node];
                    for i in 1..previous_column_nodes.len() {
                        let cur_coord = gr.node_weight(previous_column_nodes[i]).unwrap().coords;
                        // gr.node_weight_mut(previous_column_nodes[i]).unwrap().coords.v =
                        //     cur_coord.lerp(&medium_v, 0.8).v;
                        nodes.push(previous_column_nodes[i]);
                    }
                    return ([nodes.to_vec(), nodes.to_vec()], Vec::new());
                }

                for current_level in 1..dati_coords.len() {
                    built_nodes.push(Vec::new());
                    for coor in dati_coords[current_level].iter() {
                        let mut new_weight: CrochetVertex = base_weight.clone();
                        let mut delta = coor.v - new_weight.coords.v;
                        if delta <= -0.5 {
                            delta = -1.0 - delta;
                        } else if delta >= 0.5 {
                            delta = 1.0 - delta;
                        }
                        new_weight = new_weight.move_right_on_helix(delta);
                        #[cfg(debug_assertions)]
                        {
                            if new_weight.height_target().tau - new_weight.coords.tau > 2.5 {
                                eprintln!(
                                    "vecchiatau {:.3} # target {:.3}",
                                    new_weight.coords.tau,
                                    new_weight.height_target().tau
                                );
                            }
                        }

                        let new_node = gr.add_node(CrochetVertex {
                            coords: coor.clone(),
                            short_row_level: current_level as i32,
                            target_tau: Some(new_weight.height_target().tau + 1.0),
                        });
                        built_nodes[current_level].push(new_node);
                    }
                }

                // li collego
                for i in 1..built_nodes.len() {
                    connect_wales(surf, gr, &built_nodes[i - 1], &built_nodes[i]);

                    for j in 1..built_nodes[i].len() {
                        gr.add_edge(
                            built_nodes[i][j - 1],
                            built_nodes[i][j],
                            CrochetEdge {
                                link_type: EdgeType::Course,
                            },
                        );
                    }
                }

                let mut ret: [Vec<NodeIndex>; 2] = [Vec::new(), Vec::new()];

                for i in 0..built_nodes.len() {
                    ret[0].push(built_nodes[i].first().unwrap().clone());
                    ret[1].push(built_nodes[i].last().unwrap().clone());
                }

                return (ret, built_nodes.last().unwrap().to_vec());
            }

            /// Collega la patch alla colonna precedente
            fn connect_patches(
                gr: &mut Graph<CrochetVertex, CrochetEdge>,
                bordi_nodi: [&Vec<NodeIndex>; 2],
                column_middle: &Vec<TauCoords>,
            ) {
                if bordi_nodi[0][1] == bordi_nodi[1][1] {
                    // è stato fatto un join di due maglie
                    return;
                }
                for i in 1..std::cmp::min(bordi_nodi[0].len(), bordi_nodi[1].len()) {
                    // Qui collego course due a due
                    gr.add_edge(
                        bordi_nodi[0][i],
                        bordi_nodi[1][i],
                        CrochetEdge {
                            link_type: EdgeType::Course,
                        },
                    );
                }

                // Se uno dei due è più lungo di 1, collego courseoverride
                if bordi_nodi[0].len() != bordi_nodi[1].len() {
                    let min = std::cmp::min(bordi_nodi[0].len(), bordi_nodi[1].len());
                    let max = std::cmp::max(bordi_nodi[0].len(), bordi_nodi[1].len());

                    if bordi_nodi[0].len() > bordi_nodi[1].len() {
                        gr.remove_edge(
                            gr.find_edge(bordi_nodi[0][min - 1], bordi_nodi[1][min - 1])
                                .unwrap(),
                        );
                        // Fine short row
                        gr.add_edge(
                            bordi_nodi[0].last().unwrap().clone(),
                            bordi_nodi[1].last().unwrap().clone(),
                            CrochetEdge {
                                link_type: EdgeType::Course,
                            },
                        );
                        gr.add_edge(
                            bordi_nodi[0][max - 2],
                            bordi_nodi[1][min - 2],
                            CrochetEdge {
                                link_type: EdgeType::CourseNextOverrideEnd,
                            },
                        );
                    } else {
                        // Inizio short row
                        gr.add_edge(
                            bordi_nodi[0].last().unwrap().clone(),
                            bordi_nodi[1].last().unwrap().clone(),
                            CrochetEdge {
                                link_type: EdgeType::CourseNextOverrideStart,
                            },
                        );
                    }
                }
            }

            fn connect_wales(
                surf: &Box<dyn UVSurface>,
                gr: &mut Graph<CrochetVertex, CrochetEdge>,
                prev_level_nodes: &Vec<NodeIndex>,
                new_nodes: &Vec<NodeIndex>,
            ) {
                if prev_level_nodes.len() == new_nodes.len() {
                    for i in 0..prev_level_nodes.len() {
                        gr.add_edge(
                            prev_level_nodes[i],
                            new_nodes[i],
                            CrochetEdge {
                                link_type: EdgeType::Wale,
                            },
                        );
                    }
                } else if prev_level_nodes.len() == 1 && new_nodes.len() == 2 {
                    gr.add_edge(
                        prev_level_nodes[0],
                        new_nodes[0],
                        CrochetEdge {
                            link_type: EdgeType::Wale,
                        },
                    );
                    gr.add_edge(
                        prev_level_nodes[0],
                        new_nodes[1],
                        CrochetEdge {
                            link_type: EdgeType::Wale,
                        },
                    );
                } else if prev_level_nodes.len() == 2 && new_nodes.len() == 4 {
                    todo!("Ancora non supportato");
                    gr.add_edge(
                        prev_level_nodes[0],
                        new_nodes[1],
                        CrochetEdge {
                            link_type: EdgeType::Wale,
                        },
                    );
                    gr.add_edge(
                        prev_level_nodes[1],
                        new_nodes[3],
                        CrochetEdge {
                            link_type: EdgeType::Wale,
                        },
                    );
                } else if prev_level_nodes.len() == 2 && new_nodes.len() == 1 {
                    gr.update_edge(
                        prev_level_nodes[0],
                        new_nodes[0],
                        CrochetEdge {
                            link_type: EdgeType::Wale,
                        },
                    );
                    gr.update_edge(
                        prev_level_nodes[1],
                        new_nodes[0],
                        CrochetEdge {
                            link_type: EdgeType::Wale,
                        },
                    );
                } else {
                    panic!("{:?} {:?}", prev_level_nodes, new_nodes);
                }
            }

            // Un buffer che scorre e crea nodi man mano
            struct BufferCreateNodes {
                //base_nodes: [NodeIndex; 3],
                columns: [Vec<TauCoords>; 6],
                nodes: [Vec<NodeIndex>; 6],
            }

            let mut buffer: BufferCreateNodes;
            {
                let mut base_nodes: [NodeIndex; 6] = [
                    self.start_node,
                    self.start_node,
                    self.start_node,
                    self.start_node,
                    self.start_node,
                    self.start_node,
                ];
                base_nodes[1] = self.find_next_course(base_nodes[0], false).unwrap();
                base_nodes[2] = self.find_next_course(base_nodes[1], false).unwrap();
                base_nodes[3] = self.find_next_course(base_nodes[2], false).unwrap();
                base_nodes[4] = self.find_next_course(base_nodes[3], false).unwrap();
                base_nodes[5] = self.find_next_course(base_nodes[4], false).unwrap();
                let mut ultimo_nodo_alto: NodeIndex = self.start_node;
                while let Some(n) = self.find_next_course(ultimo_nodo_alto, false) {
                    ultimo_nodo_alto = n;
                }
                for _i in 0..1 {
                    let coord = self.graph.node_weight(base_nodes[0]).unwrap().coords;
                    let n = self.graph.add_node(CrochetVertex {
                        coords: coord + TauCoords { tau: 1.0, v: 0.0 },
                        short_row_level: 0,
                        target_tau: Some(2.0),
                    });
                    self.graph.add_edge(
                        base_nodes[0],
                        n,
                        CrochetEdge {
                            link_type: EdgeType::Wale,
                        },
                    );
                    self.graph.add_edge(
                        ultimo_nodo_alto,
                        n,
                        CrochetEdge {
                            link_type: EdgeType::Course,
                        },
                    );
                }
                let columns: [Vec<TauCoords>; 6] = [
                    nodes_to_generate_height(
                        self.get_surface(),
                        self.graph.node_weight(base_nodes[0]).unwrap().coords,
                        self.graph.node_weight(base_nodes[0]).unwrap().coords
                            + TauCoords { tau: 1.0, v: 0.0 },
                        false,
                    )
                    .unwrap(),
                    nodes_to_generate_height(
                        self.get_surface(),
                        self.graph.node_weight(base_nodes[1]).unwrap().coords,
                        self.graph.node_weight(base_nodes[1]).unwrap().coords
                            + TauCoords { tau: 1.0, v: 0.0 },
                        false,
                    )
                    .unwrap(),
                    nodes_to_generate_height(
                        self.get_surface(),
                        self.graph.node_weight(base_nodes[2]).unwrap().coords,
                        self.graph.node_weight(base_nodes[2]).unwrap().coords
                            + TauCoords { tau: 1.0, v: 0.0 },
                        false,
                    )
                    .unwrap(),
                    nodes_to_generate_height(
                        self.get_surface(),
                        self.graph.node_weight(base_nodes[3]).unwrap().coords,
                        self.graph.node_weight(base_nodes[3]).unwrap().coords
                            + TauCoords { tau: 1.0, v: 0.0 },
                        false,
                    )
                    .unwrap(),
                    nodes_to_generate_height(
                        self.get_surface(),
                        self.graph.node_weight(base_nodes[4]).unwrap().coords,
                        self.graph.node_weight(base_nodes[4]).unwrap().coords
                            + TauCoords { tau: 1.0, v: 0.0 },
                        false,
                    )
                    .unwrap(),
                    nodes_to_generate_height(
                        self.get_surface(),
                        self.graph.node_weight(base_nodes[5]).unwrap().coords,
                        self.graph.node_weight(base_nodes[5]).unwrap().coords
                            + TauCoords { tau: 1.0, v: 0.0 },
                        false,
                    )
                    .unwrap(),
                ];

                buffer = BufferCreateNodes {
                    columns: columns,
                    nodes: [
                        vec![
                            base_nodes[0],
                            self.find_next_course(ultimo_nodo_alto, true).unwrap(),
                        ],
                        vec![base_nodes[1]],
                        vec![base_nodes[2]],
                        vec![base_nodes[3]],
                        vec![base_nodes[4]],
                        vec![base_nodes[5]],
                    ],
                };
            }
            fn fill_future_shortrows(
                old_gr: &Graph<CrochetVertex, CrochetEdge>,
                surf: &Box<dyn UVSurface>,
                current_node: NodeIndex,
                future_shortrows: &mut Vec<(NodeIndex, usize)>,
            ) {
                const SOGLIA_ATTIVAZIONE: f64 = 3.0;
                const SOGLIA_TOLLERANZA: f64 = 1.5;

                if future_shortrows.len() > 0 {
                    return;
                }

                let coord_down = old_gr.node_weight(current_node).unwrap().coords;
                let coord_up = old_gr.node_weight(current_node).unwrap().height_target();

                let maglie: f64 = (coord_up.tau - coord_down.tau).signum()
                    * surf.geodesic_dist(coord_down, coord_up, None).unwrap()
                    / crochet_consts::get_height();

                if maglie > SOGLIA_ATTIVAZIONE {
                    let mut ultimo_nodo_pienamente_valido: usize = 0;
                    let mut ultimo_nodo_analizzato: NodeIndex = current_node;
                    // Riempio
                    loop {
                        let coord_down = old_gr.node_weight(ultimo_nodo_analizzato).unwrap().coords;
                        let coord_up = old_gr
                            .node_weight(ultimo_nodo_analizzato)
                            .unwrap()
                            .height_target();

                        let maglie: f64 = (coord_up.tau - coord_down.tau).signum()
                            * surf.geodesic_dist(coord_down, coord_up, None).unwrap()
                            / crochet_consts::get_height();

                        future_shortrows.push((ultimo_nodo_analizzato, 4));

                        if maglie > SOGLIA_ATTIVAZIONE {
                            ultimo_nodo_pienamente_valido = future_shortrows.len() - 1;
                        } else if maglie < SOGLIA_TOLLERANZA {
                            // Devo uscire
                            while future_shortrows.len() - 1 > ultimo_nodo_pienamente_valido {
                                future_shortrows.pop();
                            }
                            break;
                        }

                        if let Some(n) = static_find_next_course(old_gr, ultimo_nodo_analizzato) {
                            ultimo_nodo_analizzato = n;
                        } else {
                            future_shortrows.clear();
                            break;
                        }
                    }
                }

                if future_shortrows.len()<4
                {
                    future_shortrows.clear();
                }
            }

            // Lanciato a inizio shortrow, analizza i nodi successivi
            // indicando quanti generarne per ogni nodo base,
            //
            // Per un NodeIndex di base, quanti nodi devo generarci sopra?
            //
            // Il vettore va utilizzato solo quando non è vuoto
            let mut future_shortrows: Vec<(NodeIndex, usize)> = Vec::new();

            let mut nodes_to_smooth: Vec<NodeIndex> = Vec::new();
            const SMOOTH_WINDOW_SIZE: usize = 5;

            assert_eq!(SMOOTH_WINDOW_SIZE % 2, 1);

            let mut last_was_joined: usize = 0;
            let mut last_was_splitted: usize = 0;

            let mut need_to_close_shortrows: bool = false;

            // Reset a 0 del dithering prima di avviare la generazione
            // delle maglie
            crochet_consts::dithering_reset_tick();

            loop {
                // if buffer.columns[0][1] != buffer.columns[1][1]
                {
                    buffer.columns[3][1] = buffer.columns[3][1].lerp(
                        &buffer.columns[2][1].lerp(
                            &buffer.columns[4][1],
                            0.5
                        ),
                        0.1
                    );
                    buffer.columns[2][1] = buffer.columns[2][1].lerp(
                        &buffer.columns[1][1].lerp(
                            &buffer.columns[3][1],
                            0.5
                        ),
                        0.1
                    );
                }


                // Finché non arrivo in cima, ripeto
                #[cfg(debug_assertions)]
                {
                    // self.dot_graph(Some("tmp_graph.dot".to_string()));
                }

                const HOP_ALLOWED: i32 = 2;
                // // Vado con shortorw alte due
                for colonna in 2..6 {
                    if buffer.columns[colonna].len() as i32 % 2 != 0 {
                        buffer.columns[colonna].pop();
                    }
                    // while buffer.columns[colonna].len() as i32 > 4
                    // {
                    //     buffer.columns[colonna].pop();
                    // }
                }

                // Applico fillshortrow isteresi
                let backup_column = buffer.columns[5].to_vec();
                if future_shortrows.contains(&(buffer.nodes[5][0], 4)) {
                    while buffer.columns[5].len() < 4 {
                        buffer.columns[5].push(
                            match self.surface.height_increment(
                                buffer.columns[5].last().unwrap().clone(),
                                crochet_consts::get_height(),
                            ) {
                                Ok(e) => e,
                                _ => {
                                    buffer.columns[5] = backup_column;
                                    break;
                                }
                            },
                        );
                    }
                    while buffer.columns[5].len() > 4 {
                        buffer.columns[5].pop();
                    }

                    if buffer.nodes[5][0] == future_shortrows.last().unwrap().0 {
                        future_shortrows.clear();
                    }
                }

                // while buffer.columns[5].len() > 4 {
                //     buffer.columns[5].pop();
                // }

                // Devo togliere un elemento a quello centrale?
                {
                    let max_middle = std::cmp::max(
                        std::cmp::max(buffer.columns[1].len(), buffer.columns[2].len()),
                        std::cmp::max(buffer.columns[3].len(), buffer.columns[4].len()),
                    );
                    let min_middle = std::cmp::min(
                        std::cmp::min(buffer.columns[1].len(), buffer.columns[2].len()),
                        std::cmp::min(buffer.columns[3].len(), buffer.columns[4].len()),
                    );
                    if max_middle > buffer.columns[0].len() && max_middle > buffer.columns[5].len()
                    {
                        let target_size =
                            std::cmp::max(buffer.columns[0].len(), buffer.columns[5].len());
                        for colonna in 1..5 {
                            while buffer.columns[colonna].len() > target_size
                                && buffer.columns[colonna].len() > 2
                            {
                                for _x in 0..HOP_ALLOWED {
                                    buffer.columns[colonna].pop();
                                }
                            }
                        }
                    } else if min_middle < buffer.columns[0].len()
                        && min_middle < buffer.columns[5].len()
                    {
                        // // Evito anche singoli fine e inizio consecutivi
                        // let target_size = std::cmp::min(buffer.columns[5].len(), buffer.columns[5].len());
                        // for colonna in 1..5
                        // {
                        //     while buffer.columns[colonna].len() < target_size
                        //     {
                        //         let column_backup = buffer.columns[colonna].clone();
                        //         for _x in 0..HOP_ALLOWED
                        //         {
                        //             buffer.columns[colonna].push(
                        //                 self.surface.height_increment(
                        //                     buffer.columns[colonna].last().unwrap().clone(),
                        //                     crochet_consts::get_height()).unwrap()
                        //             );
                        //         }
                        //     }
                        // }
                    }
                }

                // Salta di più di HOP_ALLOWED alla volta?
                if buffer.columns[1].len() as i32 - buffer.columns[2].len() as i32 > HOP_ALLOWED {
                    // Se il primo è più alto del secondo
                    let diff = buffer.columns[1].len() as i32 - buffer.columns[2].len() as i32;
                    for _i in 0..diff / HOP_ALLOWED {
                        buffer.columns[1].pop();
                    }
                } else if buffer.columns[2].len() as i32 - buffer.columns[1].len() as i32
                    > HOP_ALLOWED
                {
                    let diff = buffer.columns[2].len() as i32 - buffer.columns[1].len() as i32;
                    for _i in 0..diff / HOP_ALLOWED {
                        buffer.columns[2].pop();
                    }
                }

                //println!("\n\n{:?}", buffer.columns);

                // Evito che ci siano shortrow i un nodo
                // Non molto raffinato ma funziona
                // if last_was_joined && buffer.columns[1].len()>buffer.columns[2].len()
                // {
                //     while buffer.columns[1].len() != buffer.columns[2].len()
                //     {
                //         buffer.columns[2].push(
                //             self.surface.height_increment(
                //                 buffer.columns[2].last().unwrap().clone(),
                //                 crochet_consts::get_height()).unwrap()
                //         );
                //     }
                // }

                // EFfettiva generazione dei nodi

                assert!(buffer.columns[2].len() % 2 == 0, "{:?}", buffer.columns[2]);
                let ret_coords = generate_current_step_coords(
                    &self.surface,
                    &buffer.columns[0],
                    &buffer.columns[1],
                    &buffer.columns[2],
                    last_was_joined,
                    last_was_splitted,
                );
                let coordsdata = ret_coords.0;
                last_was_joined = ret_coords.1;
                last_was_splitted = ret_coords.2;

                let (colonne, up_generated_nodes) = generate_nodes(
                    &mut self.graph,
                    &self.surface,
                    &coordsdata,
                    buffer.nodes[1][0],
                    &buffer.nodes[0],
                );
                // self.dot_graph(Some("tmp_graph.dot".to_string()));
                connect_patches(
                    &mut self.graph,
                    [&buffer.nodes[0], &colonne[0]],
                    &buffer.columns[1],
                );
                #[cfg(debug_assertions)]
                {
                    // self.dot_graph(Some("tmp_graph.dot".to_string()));
                    //self.to_string_vertices(Some("dati_grafici/oca_celeste/vertici.mat".to_string()));
                }
                buffer.nodes[1] = colonne[1].to_vec();

                // Imposto a short_row_level = 0  il nodo di base
                self.graph
                    .node_weight_mut(buffer.nodes[2][0])
                    .unwrap()
                    .short_row_level = 0;

                for n in up_generated_nodes {
                    if !nodes_to_smooth.contains(&n) {
                        nodes_to_smooth.push(n);
                    }
                }
                while nodes_to_smooth.len() > SMOOTH_WINDOW_SIZE {
                    nodes_to_smooth.remove(0);
                }

                buffer.nodes[0] = buffer.nodes[1].clone();
                buffer.nodes[1] = buffer.nodes[2].clone();
                buffer.nodes[2] = buffer.nodes[3].clone();
                buffer.nodes[3] = buffer.nodes[4].clone();
                buffer.nodes[4] = buffer.nodes[5].clone();
                buffer.nodes[5] = vec![match self.find_next_course(buffer.nodes[5][0], true) {
                    Some(el) => el,
                    None => break,
                }];
                yarn_relaxation(&mut self.graph, &self.surface, &mut nodes_to_smooth);
                // yarn_relaxation(&mut self.graph, &self.surface, &mut nodes_to_smooth);

                // All'ultimo giro devo assicurarmi di non avere short row rimaste a metà
                if !need_to_close_shortrows {
                    need_to_close_shortrows = self
                        .get_graph()
                        .node_weight(buffer.nodes[2][0])
                        .unwrap()
                        .target_tau
                        .unwrap()
                        > self.get_surface().tau_max_value() - 4.0;
                }

                // Motivo ignoto
                // un nodo è senza Wale
                // metto una pezza
                if self.find_next_wale(buffer.nodes[0][0]).is_none() {
                    self.graph.add_edge(
                        buffer.nodes[0][0],
                        buffer.nodes[0][1],
                        CrochetEdge {
                            link_type: EdgeType::Wale,
                        },
                    );
                }
                assert!(
                    self.find_next_wale(buffer.nodes[0][0]).is_some(),
                    "\n\nNo Wale a\nBuffernodes: {:.3?}\n\nBuffer columns: {:.3?}",
                    buffer.nodes,
                    buffer.columns
                );

                if !need_to_close_shortrows && crochet_consts::get_double_threshold() {
                    fill_future_shortrows(
                        self.get_graph(),
                        self.get_surface(),
                        buffer.nodes[5][0],
                        &mut future_shortrows,
                    );
                } else {
                    future_shortrows.clear();
                }

                let end_target: TauCoords = self
                    .get_graph()
                    .node_weight(buffer.nodes[2][0])
                    .unwrap()
                    .height_target();

                buffer.columns[0] = buffer.columns[1].to_vec();
                buffer.columns[1] = buffer.columns[2].to_vec();
                buffer.columns[2] = buffer.columns[3].to_vec();
                buffer.columns[3] = buffer.columns[4].to_vec();
                buffer.columns[4] = buffer.columns[5].to_vec();
                buffer.columns[5] = match nodes_to_generate_height(
                    self.get_surface(),
                    self.get_graph()
                        .node_weight(buffer.nodes[5][0])
                        .unwrap()
                        .coords,
                    end_target,
                    need_to_close_shortrows,
                ) {
                    Ok(e) => e,
                    Err(_) => break,
                };
                if need_to_close_shortrows && buffer.columns[2].len() > 2
                    || !crochet_consts::generate_short_row_flag()
                {
                    let new_column_size = match buffer.columns[1].len() > 2 {
                        true => buffer.columns[1].len() - 1,
                        false => 2,
                    };
                    while buffer.columns[2].len() > new_column_size {
                        buffer.columns[2].pop();
                    }
                }
                //assert!(end_target.tau - buffer.columns[2][0].tau > 0.05);
                //test_check_wale_dist(&self.surface, &buffer.columns[2]);

                buffer.columns[2][0] = self.graph.node_weight(buffer.nodes[2][0]).unwrap().coords;
                #[cfg(debug_assertions)]
                {
                    // println!("passaggio prossimo nodo: {:.3?}", buffer.columns[2][0]);
                    // self.to_string_vertices(Some("dati_grafici/godot_meshviz/vertici_gismondi.mat".to_string()));
                    // self.debug_shortrow_level_nodes(Some("dati_grafici/godot_meshviz/vertici_gismondi_shortrow1.mat".to_string()), 1);
                }

                {
                    // Adatto minimo tra tau e target_tau a 1
                    let w_0: &CrochetVertex =
                        self.get_graph().node_weight(buffer.nodes[0][0]).unwrap();
                    let w_1: &CrochetVertex =
                        self.get_graph().node_weight(buffer.nodes[1][0]).unwrap();
                    if w_0.coords.v < 0.5 && w_1.coords.v >= 0.5
                    ||
                    w_0.coords.v < 0.0 && w_1.coords.v >= 0.0
                    {
                        // Una volta a giro. Intorno a 0.5, v è monotona strettamente crescente
                        // Trovo il minimo di delta rispetto al target
                        let mut reached_top_level: bool = false;
                        let mut min_found: f64 = f64::MAX;
                        let mut start_ind = buffer.nodes[0][0];
                        let mut fuond_we = self.get_graph().node_weight(start_ind).unwrap();
                        while let Some(new_ind) = self.find_next_course(start_ind, true) {
                            let w: &CrochetVertex =
                                self.get_graph().node_weight(start_ind).unwrap();
                            let maglie = match self.get_surface().geodesic_dist(
                                w.coords,
                                w.height_target(),
                                None,
                            ) {
                                Ok(e) => e,
                                _ => {
                                    reached_top_level = true;
                                    break;
                                }
                            } / crochet_consts::get_height();
                            if maglie < min_found {
                                min_found = maglie;
                                fuond_we = w;
                            }

                            start_ind = new_ind;
                        }

                        // Applico correzione
                        let correzione_tau: f64 = if min_found > 1.8 { 1.0 } else { 0.0 };
                        let mut start_ind = buffer.nodes[0][0];
                        while let Some(new_ind) = self.find_next_course(start_ind, true) {
                            if reached_top_level {
                                break;
                            }
                            let w: &mut CrochetVertex =
                                self.graph.node_weight_mut(start_ind).unwrap();
                            w.target_tau = Some(w.target_tau.unwrap() - correzione_tau);

                            start_ind = new_ind;
                        }
                        let w: &mut CrochetVertex = self.graph.node_weight_mut(start_ind).unwrap();
                        w.target_tau = Some(w.target_tau.unwrap() - correzione_tau);
                    }
                }

                // Avanzo il tick ciclicamente
                // una volta per iterazione
                crochet_consts::dithering_advance_tick();
            }

            {
                // Serve a smussare gli ultimi nodi aggiunti in alto alla fine di tutto
                // E impostarli a 0 di short_row_level
                let mut last_node = buffer.nodes[1][0];
                while let Some(tmp) = self.find_next_course(last_node, true) {
                    self.graph.node_weight_mut(tmp).unwrap().short_row_level = 0;

                    // yarn_relaxation(&mut self.graph, &self.surface, last_node);

                    last_node = tmp;
                }
            }
            self.fix_graph();
            self.undeform();
            // self.undeform();
            self.test_graph_validity();
        }

        fn get_start_node(&self) -> NodeIndex {
            return self.start_node.clone();
        }
        fn get_graph<'a>(&'a self) -> &'a Graph<CrochetVertex, CrochetEdge> {
            &self.graph
        }
        fn get_surface<'a>(&'a self) -> &'a Box<dyn UVSurface> {
            &self.surface
        }
    }
}
