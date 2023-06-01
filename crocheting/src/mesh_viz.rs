pub mod mesh_viz {
    use crate::nodes::nodes::{EdgeType, KnittableGraph};
    use crate::param_surface::surface::{UVSurface, TauCoords, Cylinder};
    use nalgebra::{Point3, Vector3};
    use petgraph::prelude::{EdgeIndex, NodeIndex, StableGraph, Direction};
    use petgraph::visit::EdgeRef;
    use crate::geometry::curves::Bezier;
    use crate::geometry::transform::TriLinearInterpolation;
    use crate::crochet_consts::crochet_consts::get_yarn_radius;
    use std::fs::File;
    use std::fs;
    use std::io::Write;
    use crate::crochet_consts;

    

    #[derive(PartialEq, Eq, Hash, Debug, Clone, Copy)]
    pub enum StitchType {
        SlipSitch = 0,     //Bassissima
        SingleCrochet = 1, //Bassa
        TurnBack = 2,      //Inverti direzione, vai indietro
        TurnForward = 3,   //Andavo indietro e torno dritto
        Increase = 4,      //Aumento
        Decrease = 5,      //Diminuzione
        EndingCap = 6      // Serve visivamente ma non deve essere usato nelle istruzioni
    }

    #[derive(Debug, PartialEq, Clone)]
    pub enum StitchConnection {
        NextWale,
        NextCourse,
        PrevWale,
        PrevCourse,
    }

    #[derive(Debug)]
    struct StitchEdge {
        link: StitchConnection,
    }

    #[derive(Debug)]
    pub struct StitchMesh {
        graph: StableGraph<Box<dyn StitchFace>, StitchEdge>,
        cache_vertices: Option<(Vec<Point3<f64>>, Vec<i32>)>,
        first_stitch: Option<NodeIndex>
    }

    impl StitchMesh {
        pub fn new() -> Self
        {
            StitchMesh {
                graph: StableGraph::new(),
                cache_vertices: None,
                first_stitch: None
            }
        }

        pub fn new_from_graph(old_gr: &Box<dyn KnittableGraph>) -> Self {
            let mut st = StitchMesh {
                graph: StableGraph::new(),
                cache_vertices: None,
                first_stitch: None
            };

            let shortrows = old_gr.shortrows_limits();

            fn next_edge(old_gr: &Box<dyn KnittableGraph>, edge: EdgeIndex, find_start: bool, find_end: bool) -> Option<EdgeIndex>
            {
                let gr = old_gr.get_graph();
                let target: NodeIndex = gr.edge_endpoints(edge).unwrap().1;
                return old_gr.find_next_course_edge(target, find_start, find_end);
            }

            /// Chiamato per ogni Edge, mi dice a quali Edge sopra è collegato
            fn up_edges(old_gr: &Box<dyn KnittableGraph>, edge: EdgeIndex) -> Result<Vec<EdgeIndex>, &str>
            {
                let gr = old_gr.get_graph();
                let endpoints = gr.edge_endpoints(edge).unwrap();

                let first_up = match old_gr.ordered_wales(endpoints.0).first()
                {
                    Some(e)=>e.clone(),
                    None=>return Err("no wale")
                };
                let end_up = match old_gr.ordered_wales(endpoints.1).last()
                {
                    Some(e)=>e.clone(),
                    None=>return Err("no wale")
                };

                if first_up == end_up
                {
                    return Ok(Vec::new());
                }

                if let Some(ed) = gr.find_edge(first_up, end_up)
                {
                    return Ok(vec![
                        ed
                    ]);
                }
                else
                {
                    #[cfg(debug_assertions)]
                    {
                        println!("nodo {}", old_gr.print_node(first_up));
                    }
                    
                    let n: EdgeIndex = match old_gr.find_next_course_edge(first_up, false, false)
                    {
                        Some(e) => e,
                        None => {return Err("No prossimo nodo a up_edge")}
                    }; 
                    let next: EdgeIndex = match next_edge(old_gr, n, false, false)
                    {
                        Some(e) => e,
                        None => {return Err("No prossimo nodo a up_edge")}
                    }; 
                    #[cfg(debug_assertions)]
                    {
                        println!("nodo {}", old_gr.print_node(first_up));
                    }
                    #[cfg(debug_assertions)]
                    {
                        println!("nodo {}", old_gr.print_node(end_up));
                    }

                    if old_gr.find_all_wales(endpoints.0).len() == 1
                    {
                        return Ok(
                            vec![
                                n
                            ]
                        )
                    }


                    return Ok(vec![
                        n,
                        next
                    ]);
                }
            }

            /// Vettore di StitchFace, già invertite di direzione
            /// 
            /// Restituisce anche il nuovo EdgeIndex da cui continuare ad analizzare
            /// 
            /// L'ultimo terzo valore della tupla indica se una faccia è l'inizio di un giro
            /// - restituito come usize → indice di quale è la faccia da usare come inizio riga
            /// - restituito come EdgeIndex, il edge buffer sotto in corrispondenza di quella faccia genereata
            fn resolve_backward_shortrow(
                old_gr: &Box<dyn KnittableGraph>,
                shortrowlimits: &Vec<(NodeIndex, NodeIndex)>,
                end_vertex: NodeIndex,
                current_rowstartedge: EdgeIndex
            ) -> (Vec<(Box<dyn StitchFace>, EdgeIndex)>, EdgeIndex, Option<(usize, EdgeIndex)>)
            {
                let mut accumulated_faces: Vec<(Box<dyn StitchFace>, EdgeIndex)> = Vec::new();

                let mut start_of_row: Option<(usize, EdgeIndex)> = None;

                let start_vertex: NodeIndex = shortrowlimits.iter().filter(
                    |e| {e.1 == end_vertex}
                )
                .map(|e|e.0)
                .collect::<Vec<NodeIndex>>()[0];

                {
                    let up_1: NodeIndex = old_gr.ordered_wales(start_vertex)[0];
                    let up_2: NodeIndex = old_gr.ordered_wales(up_1)[0];
                    let root_overr: NodeIndex = old_gr.get_graph()
                    .edges_directed(start_vertex, Direction::Incoming)
                    .filter(
                        |e|{e.weight().link_type == EdgeType::Course}
                    ).next().unwrap().source();
                    // aggiungo il primissimo turnforward


                    {
                        // Calcolo se è un inizio giro
                        let latosotto = old_gr.get_graph().find_edge(
                            root_overr,
                            start_vertex
                        ).unwrap();
                        if current_rowstartedge == latosotto
                        {
                            start_of_row = Some(
                                (
                                    0,
                                    old_gr.get_graph().find_edge(
                                        root_overr,
                                        up_2
                                    ).unwrap()
                                )
                            )
                        }
                    }



                    let node_indices = [
                        root_overr,
                        start_vertex,
                        up_1,
                        up_2
                    ];


                    let node_coords: [TauCoords; 4] = node_indices.iter().map(
                        |e|{old_gr.node_taucoords(*e)}
                    ).collect::<Vec<TauCoords>>().try_into().unwrap();

                    accumulated_faces.push(
                        (
                            TurnForwardFace::new(
                                node_coords,
                                old_gr.get_surface()
                            ),
                            old_gr.get_graph().find_edge(root_overr, up_2).unwrap()
                        )
                        
                    );
                }
                

                // Analizzo e scorro tra start_vertex ed end_vertex
                let mut moving_analyzer: Vec<EdgeIndex>;
                if let Some(n2) = old_gr.filter_edges(start_vertex, EdgeType::Course).first()
                {
                    let ed1 = old_gr.get_graph().find_edge(start_vertex, n2.clone()).unwrap();
                    moving_analyzer = vec![
                        ed1
                    ];
                    if let Some(e2) = next_edge(old_gr, ed1, false, false)
                    {
                        moving_analyzer.push(e2);
                    }
                }
                else
                {
                    panic!("Finestra shortrow non inizializzabile")
                }

                loop
                {
                    let mut double_skip = false;
                    if let Ok(res_face) = create_stitchface(old_gr, &moving_analyzer)
                    {
                        let face = &res_face.0;
                        // Se invece è una faccia normale, la aggiungo
                        if face.get_stitch_type() == StitchType::Decrease
                        {
                            //Se questa faccia normale è una decrease, attivo double_slip
                            double_skip = true;
                        }


                        {
                            // Calcolo se è un inizio giro
                            if current_rowstartedge == moving_analyzer[0]
                            ||
                            (
                                double_skip
                                &&
                                current_rowstartedge == moving_analyzer[1]
                            )
                            {
                                start_of_row = Some(
                                    (
                                        accumulated_faces.len(),
                                        res_face.1
                                    )
                                )
                            }
                        }


                        accumulated_faces.push(res_face);
                        
                    }
                    else
                    {
                        break;
                    }


                    // Avanzo finestra
                    if let Some(e) = next_edge(old_gr, moving_analyzer.last().unwrap().clone(), false, false)
                    {
                        moving_analyzer.push(e);
                    }
                    moving_analyzer.remove(0);
                    if moving_analyzer.len() == 0
                    {
                        break;
                    }
                
                    if double_skip
                    {
                        
                        if let Some(e) = next_edge(old_gr, moving_analyzer.last().unwrap().clone(), true, true)
                        {
                            moving_analyzer.push(e);
                        }
                        moving_analyzer.remove(0);
                        if moving_analyzer.len() == 0
                        {
                            break;
                        }
                    }
                }

                accumulated_faces.reverse();
                if let Some(s) = start_of_row
                {
                    start_of_row = Some((accumulated_faces.len()-1 - s.0, s.1));
                }
                
                return
                (
                    accumulated_faces,
                    old_gr.find_next_course_edge(
                        old_gr.ordered_wales(start_vertex)[0],
                        false, false
                    ).unwrap(),
                    start_of_row
                );
            }

            /// Prende i lati attualmente nei buffer e crea la stitchface corretta da aggiungere
            /// 
            /// Restituisce anche l'EdgeIndex sovrastante, da usare per rilevare l'inizio giro
            fn create_stitchface(old_gr: &Box<dyn KnittableGraph>, edge_low_buffer: &Vec<EdgeIndex>) -> Result<(Box<dyn StitchFace>, EdgeIndex), String>
            {
                /// Ci troviamo in una fine di una shortrow?
                fn is_end_turn(old_gr: &Box<dyn KnittableGraph>, edge: EdgeIndex) -> bool
                {
                    return
                    old_gr.get_graph().edges_directed(
                        old_gr.get_graph().edge_endpoints(edge).unwrap().1,
                        Direction::Incoming
                        )
                        .filter(|e|{e.weight().link_type == EdgeType::CourseNextOverrideEnd})
                        .count() == 1
                        &&
                        old_gr.get_graph().edge_weight(edge).unwrap().link_type != EdgeType::CourseNextOverrideEnd;
                }
                fn is_increase( edges_up: &Vec<Vec<EdgeIndex>>) -> bool
                {
                    return edges_up[0].len() > 1;
                }
                fn is_decrease(edges_low: &Vec<EdgeIndex>, edges_up: &Vec<Vec<EdgeIndex>>) -> bool
                {
                    if edges_low.len() != 2
                    {
                        return false;
                    }
                    if edges_up.len()>1 && edges_up[0].len()==0 && edges_up[1].len()==1
                    {
                        return true;
                    }
                    else
                    {
                        return false;
                    }
                }
                fn is_singlecrochet(edges_low: &Vec<EdgeIndex>, edges_up: &Vec<Vec<EdgeIndex>>) -> bool
                {
                    if edges_low.len()==1
                    {
                        if edges_up[0].len()==1
                        {
                            return true;
                        }
                    }
                    else
                    {
                        if edges_up[0].len()!=1
                        {
                            // Sicuro un aumento
                            return false
                        }
                        // Il grosso dovrebbe essere stato filtrato prima
                        // dalle altre funzioni chiamate
                        return true;
                    }
                    return false;
                }
                
                if is_end_turn(old_gr, edge_low_buffer[0])
                {
                    // Fine shortrow

                    // Faccio prima di tutto l'end turn perché altrimenti
                    // mentre calcola gli up_edges è facile che esplode perché manca un Course che si aspetta

                    let end_points_low = old_gr.get_graph().edge_endpoints(edge_low_buffer[0]).unwrap();

                    // DEVO PRENDERE I WALE IN ORDINE
                    let up_1 = old_gr.ordered_wales(end_points_low.0).last().unwrap().clone();
                    let up_2 = old_gr.ordered_wales(up_1).last().unwrap().clone();

                    let node_indices = [
                        end_points_low.0,
                        end_points_low.1,
                        up_2,
                        up_1
                    ];
                    let node_coords: [TauCoords; 4] = node_indices.iter().map(
                        |e|{old_gr.node_taucoords(*e)}
                    ).collect::<Vec<TauCoords>>().try_into().unwrap();

                    return Ok(
                        (
                            TurnBackFace::new(
                                node_coords,
                                old_gr.get_surface()
                            ),
                            old_gr.get_graph().find_edge(up_2, end_points_low.1).unwrap()
                        ));
                }

                let comp_up_edges: Vec<Vec<EdgeIndex>>;
                {
                    let mut accumul: Vec<Vec<EdgeIndex>> = Vec::new();
                    let compup_with_errors: Vec<Result<Vec<EdgeIndex>, &str>> = edge_low_buffer.iter().map(|e|{up_edges(old_gr, *e)}).collect();
                    for (i, el) in compup_with_errors.iter().enumerate()
                    {
                        if el.is_err()
                        {
                            if i == 0
                            {
                                return Err(el.clone().unwrap_err().into())
                            }
                            
                        }
                        else
                        {
                            accumul.push(el.clone().unwrap());
                        }
                    }
                    comp_up_edges = accumul;
                }
                if is_increase(&comp_up_edges)
                {
                    // Aumento
                    let end_points_1_low = old_gr.get_graph().edge_endpoints(edge_low_buffer[0]).unwrap();

                    let up_endpoints_1 = old_gr.get_graph().edge_endpoints(comp_up_edges[0][0]).unwrap();
                    let up_endpoints_2 = old_gr.get_graph().edge_endpoints(comp_up_edges[0][1]).unwrap();

                    let node_indices = [
                        end_points_1_low.0,
                        end_points_1_low.1,
                        up_endpoints_2.1,
                        up_endpoints_2.0,
                        up_endpoints_1.0
                    ];
                    let node_coords: [TauCoords; 5] = node_indices.iter().map(
                        |e|{old_gr.node_taucoords(*e)}
                    ).collect::<Vec<TauCoords>>().try_into().unwrap();

                    Ok(
                        (
                            IncreaseFace::new(
                                node_coords,
                                old_gr.get_surface()
                            ),
                            comp_up_edges[0][1]
                        )
                    )
                }
                else if is_decrease(edge_low_buffer, &comp_up_edges)
                {
                    // Diminuzione
                    let end_points_1_low = old_gr.get_graph().edge_endpoints(edge_low_buffer[0]).unwrap();
                    let end_points_2_low = old_gr.get_graph().edge_endpoints(edge_low_buffer[1]).unwrap();

                    let up_endpoints = old_gr.get_graph().edge_endpoints(comp_up_edges[1][0]).unwrap();

                    let node_indices = [
                        end_points_1_low.0,
                        end_points_1_low.1,
                        end_points_2_low.1,
                        up_endpoints.1,
                        up_endpoints.0
                    ];
                    let node_coords: [TauCoords; 5] = node_indices.iter().map(
                        |e|{old_gr.node_taucoords(*e)}
                    ).collect::<Vec<TauCoords>>().try_into().unwrap();
                    Ok(
                        (
                            DecreaseFace::new(
                                node_coords,
                                old_gr.get_surface()
                            ),
                            comp_up_edges[1][0]
                        )
                    )


                    
                }
                else if is_singlecrochet(edge_low_buffer, &comp_up_edges)
                {
                    // Maglia bassa
                    let end_points_low = old_gr.get_graph().edge_endpoints(edge_low_buffer[0]).unwrap();

                    let up_endpoints = old_gr.get_graph().edge_endpoints(comp_up_edges[0][0]).unwrap();

                    let node_indices = [
                        end_points_low.0,
                        end_points_low.1,
                        up_endpoints.1,
                        up_endpoints.0
                    ];
                    let node_coords: [TauCoords; 4] = node_indices.iter().map(
                        |e|{old_gr.node_taucoords(*e)}
                    ).collect::<Vec<TauCoords>>().try_into().unwrap();

                    Ok(
                        (
                            SingleCrochetFace::new(
                                node_coords,
                                old_gr.get_surface()
                            ),
                            comp_up_edges[0][0]
                        )
                    )
                }
                else
                {
                    return Err("non trovato".into());
                }
            
            }


            let mut moving_analyzer: Vec<EdgeIndex>;
            if let Some(n2) = old_gr.filter_edges(old_gr.get_start_node(), EdgeType::Course).first()
            {
                let ed1 = old_gr.get_graph().find_edge(old_gr.get_start_node(), n2.clone()).unwrap();
                moving_analyzer = vec![
                    ed1,
                    next_edge(old_gr, ed1, false, false).unwrap()
                ];
            }
            else
            {
                return st;
            }


            let mut last_face_added: Option<NodeIndex> = None;
            
            // Mi serve per costruire il filone di wale per 
            // la separazione dei giri
            // come
            // (face_index, vecchiograph_edge)
            let mut first_of_row: Option<(NodeIndex, EdgeIndex)> = None;

            // Sposto moving_analyzer
            // e creo tutti i nodi "faccia"
            loop
            {
                let mut double_skip = false;
                if let Ok(res_face) = create_stitchface(old_gr, &moving_analyzer)
                {
                    let face = res_face.0;
                    let main_up_edge: EdgeIndex = res_face.1;
                    

                    // Se è una fine shortrow
                    if face.get_stitch_type() == StitchType::TurnBack
                    {
                        let new_face_index = st.add(face, last_face_added);
                        last_face_added = Some(new_face_index);

                        if st.connect_wale(first_of_row, moving_analyzer[0], new_face_index)
                        {
                            first_of_row = Some(
                                (
                                    new_face_index,
                                    main_up_edge
                                )
                            );
                        }
                        
                        // lancio calcolo shortrow e le aggiungo
                        let res = resolve_backward_shortrow(old_gr, &shortrows,
                            old_gr.ordered_wales(
                                old_gr.get_graph().edge_endpoints(moving_analyzer[0]).unwrap().0
                            ).last().unwrap().clone(),
                            first_of_row.unwrap().1
                        );
                        
                        let mut k: usize = 0;
                        for f in res.0
                        {
                            let new_face_index = st.add(f.0, last_face_added);                            
                            last_face_added = Some(new_face_index);

                            if let Some(rowstart_added) = res.2
                            {
                                if rowstart_added.0 == k
                                {
                                    st.graph.add_edge(first_of_row.unwrap().0, new_face_index, StitchEdge{link: StitchConnection::NextWale});
                                    first_of_row = Some(
                                        (
                                            new_face_index,
                                            rowstart_added.1
                                        )
                                    );
                                }
                            }
                            k += 1;
                        }
                        

                        // Cambio i dati nel buffer lati in basso
                        moving_analyzer = vec![res.1];
                        if let Some(e) = next_edge(old_gr, moving_analyzer.last().unwrap().clone(), true, true)
                        {
                            moving_analyzer.push(e);
                        }
                        continue;
                    }
                    else
                    {
                        // Se invece è una faccia normale, la aggiungo
                        let tipo = face.get_stitch_type();
                        if tipo == StitchType::Decrease
                        {
                            //Se questa faccia normale è una decrease, attivo double_slip
                            double_skip = true;
                        }
                        
                        let new_face_index = st.add(face, last_face_added);
                        last_face_added = Some(new_face_index);
                        
                        if st.connect_wale(first_of_row, moving_analyzer[0], new_face_index)
                            ||
                            (
                                tipo == StitchType::Decrease
                                &&
                                st.connect_wale(first_of_row, moving_analyzer[1], new_face_index)
                            )
                        {
                            first_of_row = Some(
                                (
                                    new_face_index,
                                    main_up_edge
                                )
                            );
                        }
                    }
                }
                else
                {
                    if old_gr.find_all_wales(old_gr.get_graph().edge_endpoints(moving_analyzer[0]).unwrap().0).len() > 0
                    {
                        // Continuo e ignoro
                    }
                    else
                    {
                        break;
                    }
                }


                // Avanzo finestra
                moving_analyzer.remove(0);
                if let Some(e) = next_edge(old_gr, moving_analyzer.last().unwrap().clone(), true, true)
                {
                    moving_analyzer.push(e);
                }
                if moving_analyzer.len() == 0
                {
                    break;
                }
                
                if double_skip
                {
                    moving_analyzer.remove(0);
                    if let Some(e) = next_edge(old_gr, moving_analyzer.last().unwrap().clone(), true, true)
                    {
                        moving_analyzer.push(e);
                    }
                    if moving_analyzer.len() == 0
                    {
                        break;
                    }
                }
            }

            
            
            
            // Aggiungo terminazinoe alta
            let mut cur_edge = moving_analyzer[0];
            // #[cfg(debug_assertions)]
            {
                println!("{}",old_gr.print_node(
                    old_gr.get_graph().edge_endpoints(cur_edge).unwrap().0
                ));
                println!("{}",old_gr.print_node(
                    old_gr.get_graph().edge_endpoints(cur_edge).unwrap().1
                ));
                if let Ok(res_face) = create_stitchface(old_gr, &moving_analyzer)
                {
                    println!("{:?}", res_face);
                }
            }
            loop {
                let cur_endpoints = old_gr.get_graph().edge_endpoints(cur_edge).unwrap();
                
                let low_low_indices: (NodeIndex, NodeIndex);
                
                if let (Some(a), Some(b)) = (
                    old_gr.get_graph().edges_directed(cur_endpoints.0, Direction::Incoming
                    ).filter(|e|{e.weight().link_type == EdgeType::Wale})
                    .map(|e|e.source())
                    .next(),
                    old_gr.get_graph().edges_directed(cur_endpoints.0, Direction::Incoming
                    ).filter(|e|{e.weight().link_type == EdgeType::Wale})
                    .map(|e|e.source())
                    .next()
                )
                {
                    low_low_indices = (a, b);
                }
                else
                {
                    break;
                }

                let surf = old_gr.get_surface();

                let cur_points: (Point3<f64>, Point3<f64>) =
                    (
                    surf.to_3d_point(
                        old_gr.node_taucoords(cur_endpoints.0)
                    ).unwrap(),
                    surf.to_3d_point(
                        old_gr.node_taucoords(cur_endpoints.1)
                    ).unwrap()
                );
                let dir: (Vector3<f64>, Vector3<f64>) = 
                (
                    cur_points.0
                    -
                    surf.to_3d_point(
                        old_gr.node_taucoords(low_low_indices.0)
                    ).unwrap(),

                    cur_points.1
                    -
                    surf.to_3d_point(
                        old_gr.node_taucoords(low_low_indices.1)
                    ).unwrap()
                );
                
                let up_points = (
                    cur_points.0
                    +
                    dir.0.normalize()*crochet_consts::crochet_consts::get_height(),
                    cur_points.1
                    +
                    dir.1.normalize()*crochet_consts::crochet_consts::get_height(),
                );

                let normals = (
                    surf.normal_versor(old_gr.node_taucoords(cur_endpoints.0)),
                    surf.normal_versor(old_gr.node_taucoords(cur_endpoints.1))
                );

                let spessore = crochet_consts::crochet_consts::get_thickness();

                let below_face = [
                    cur_points.0 - normals.0*spessore,
                    cur_points.1 - normals.1*spessore,
                    up_points.1 - normals.1*spessore,
                    up_points.0 - normals.0*spessore,
                ];
                let up_face = [
                    cur_points.0 + normals.0*spessore,
                    cur_points.1 + normals.1*spessore,
                    up_points.1 + normals.1*spessore,
                    up_points.0 + normals.0*spessore,
                ];

                
                st.add(CapFace::new(
                    below_face,
                    up_face
                ), None);

                if let Some(e) = next_edge(old_gr, cur_edge, true, true)
                {
                    cur_edge = e;
                }
                else
                {
                    break;
                }
            }


            return st;
        }

        fn add(&mut self, f: Box<dyn StitchFace>, prev_face: Option<NodeIndex>) -> NodeIndex
        {
            let new_face = self.graph.add_node(f);
            if let Some(last_f) = prev_face
            {
                self.connect_course(last_f, new_face);
            }
            if self.first_stitch.is_none()
            {
                self.first_stitch = Some(new_face);
            }
            return new_face
        }

        fn connect_course(&mut self, f_prev: NodeIndex, f_next: NodeIndex) -> ()
        {
            self.graph.add_edge(f_prev, f_next, StitchEdge{link: StitchConnection::NextCourse});
        }

        /// Se necessario, aggiunge un wale
        /// 
        /// Restituisce i nuovi sì o no, a seconda se qualcosa è stato aggiunto
        fn connect_wale(
            &mut self,
            last_startrow_data: Option<(NodeIndex, EdgeIndex)>,
            cur_buffer_edge: EdgeIndex,
            new_face: NodeIndex
        ) -> bool
        {
            if last_startrow_data.is_none()
            {
                return true;
            }




            if let Some(data) = last_startrow_data
            {
                if data.1 == cur_buffer_edge
                {
                    self.graph.add_edge(data.0, new_face, StitchEdge{link: StitchConnection::NextWale});
                    return true;
                }
            }

            return false;
        }

        /// Salva su file vertici e indici delle facce senza triangolarle.
        /// 
        /// Il formato risultante è un elenco delle coordinate dei vertici di ogni faccia
        /// Ignoro la gestione degli indici, che saranno semplicemente incrementati uno a uno
        pub fn data_for_blender_viz(&self, filename: &str)
        {
            let mut output = File::create(filename).unwrap();
            for ind in self.graph.node_indices()
            {   
                let boundaries: [Vec<Point3<f64>>;2] = self.graph.node_weight(ind).unwrap().boundary_vertices();
                for bound in boundaries
                {
                    for (j, point3) in bound.iter().enumerate()
                    {
                        write!(output, "{} {} {}", point3.x, point3.y, point3.z).unwrap();
    
                        if j != bound.len()-1
                        {
                            write!(output, "|").unwrap();
                        }
                    }
                    if bound.len()!=0
                    {
                        write!(output, "\n").unwrap();
                    }
                }
                
            }
        }

        ///
        /// Il parametro modifica la visualizzazione da "mesh semplice struttura" a "mesh completa delle maglie"
        fn generate_vertex_cache(&mut self, complete_stitches_visual: bool)
        {
            if self.cache_vertices.is_none()
            {
                let mut end_result: (Vec<Point3<f64>>, Vec<i32>) = (Vec::new(), Vec::new());
                let mut starting_index = 0;
                for ind in self.graph.node_indices()
                {
                    let res_array: (Vec<Point3<f64>>, Vec<i32>);
                    
                    if complete_stitches_visual
                    {
                        res_array = self.graph.node_weight(ind).unwrap().all_vertex_data(starting_index);
                    }
                    else
                    {
                        res_array = self.graph.node_weight(ind).unwrap().mesh_only_vertex_data(starting_index);
                    }
                    

                    starting_index = *res_array.1.iter().max().unwrap() + 1;

                    end_result.0.extend(res_array.0);
                    end_result.1.extend(res_array.1);
                }
                self.cache_vertices = Some(
                    end_result
                )
            }
        }

        pub fn array_verts(&mut self, complete_stitches_visual: bool) -> Vec<Point3<f64>>
        {
            self.generate_vertex_cache(complete_stitches_visual);
            return self.cache_vertices.as_ref().unwrap().0.clone();
        }
        pub fn array_vertex_indices(&mut self, complete_stitches_visual: bool) -> Vec<i32>
        {
            self.generate_vertex_cache(complete_stitches_visual);
            return self.cache_vertices.as_ref().unwrap().1.clone();
        }

        pub fn all_face_data(&self) -> Vec<(StitchType, Vec<Point3<f64>>, Vec<Point3<f64>>)>
        {
            let mut v: Vec<(StitchType, Vec<Point3<f64>>, Vec<Point3<f64>>)> = Vec::new();

            for el in self.graph.node_indices() {
                v.extend(self.graph.node_weight(el).unwrap().face_data());
            }
            return v;
        }

        pub fn all_face_data_to_file(&self, filename: &str)
        {
            let dati = self.all_face_data();
            
            let mut output = File::create(filename).unwrap();
            for faccia in dati
            {
                write!(output, "{:?}\n", faccia.0).unwrap();

                for v in faccia.1
                {
                    write!(output, "{} {} {}\n", v.x, v.y, v.z).unwrap();
                }
                for v in faccia.2
                {
                    write!(output, "{} {} {}\n", v.x, v.y, v.z).unwrap();
                }
                write!(output, "\n").unwrap();
            }
        }

        pub fn get_trace(&self, connection: StitchConnection) -> Vec<(NodeIndex, StitchType)>
        {
            fn next(graph: &StableGraph<Box<dyn StitchFace>, StitchEdge>, start: NodeIndex, conn: StitchConnection) -> Option<(NodeIndex, StitchType)>
            {
                let b: Vec<(NodeIndex, StitchType)> = graph
                .edges_directed(start, Direction::Outgoing)
                .filter(|e|{e.weight().link == conn})
                .map(|e|{e.target()})
                .map(|i|{(i, graph.node_weight(i).unwrap().get_stitch_type())})
                .collect();
                if b.len()==0
                {
                    return None
                }
                else if b.len()==1
                {
                    return Some(b[0]);
                }
                else
                {
                    panic!("Due course stitch")
                }
            }

            let mut faces: Vec<(NodeIndex, StitchType)> = vec![
                (self.first_stitch.unwrap(), self.graph.node_weight(self.first_stitch.unwrap()).unwrap().get_stitch_type() )
                ];
            while let Some(newnode) = next(&self.graph, faces.last().unwrap().0.clone(), connection.clone())
            {
                faces.push(newnode);
            }

            return faces;
        }

        pub fn node_count(&self) -> usize
        {
            return self.graph.node_count();
        }
    }
  


    // Curve
    lazy_static!
    {
        static ref SINGLE_STITCH_CURVES: Vec<Vec<(Point3<f64>, Point3<f64>,Point3<f64>)>> = 
        read_curve_data("SingleCrochet_curvedata.txt");
    }
    lazy_static!
    {
        static ref INCREASE_L_STITCH_CURVES: Vec<Vec<(Point3<f64>, Point3<f64>,Point3<f64>)>> = 
        read_curve_data("IncreaseLeft_curvedata.txt");
    }
    lazy_static!
    {
        static ref INCREASE_R_STITCH_CURVES: Vec<Vec<(Point3<f64>, Point3<f64>,Point3<f64>)>> = 
        read_curve_data("IncreaseRight_curvedata.txt");
    }
    lazy_static!
    {
        static ref DECREASE_L_STITCH_CURVES: Vec<Vec<(Point3<f64>, Point3<f64>,Point3<f64>)>> = 
        read_curve_data("DecreaseLeft_curvedata.txt");
    }
    lazy_static!
    {
        static ref DECREASE_R_STITCH_CURVES: Vec<Vec<(Point3<f64>, Point3<f64>,Point3<f64>)>> = 
        read_curve_data("DecreaseRight_curvedata.txt");
    }
    lazy_static!
    {
        static ref TURNBACK_STITCH_CURVES: Vec<Vec<(Point3<f64>, Point3<f64>,Point3<f64>)>> = 
        read_curve_data("TurnBack_curvedata.txt");
    }
    lazy_static!
    {
        static ref TURNFORWARD_STITCH_CURVES: Vec<Vec<(Point3<f64>, Point3<f64>,Point3<f64>)>> = 
        read_curve_data("TurnForward_curvedata.txt");
    }
    lazy_static!
    {
        static ref CAP_STITCH_CURVES: Vec<Vec<(Point3<f64>, Point3<f64>,Point3<f64>)>> = 
        read_curve_data("Cap_curvedata.txt");
    }


    fn read_curve_data(filename: &str) -> Vec<Vec<(Point3<f64>, Point3<f64>,Point3<f64>)>>
    {
        let mut contenuto: String = fs::read_to_string(filename).unwrap();
        contenuto = contenuto.trim_end().to_string();

        let mut ret: Vec<Vec<(Point3<f64>, Point3<f64>,Point3<f64>)>> = Vec::new();
        for (i, d) in contenuto.split("\n\n").enumerate()
        {
            ret.push(Vec::new());
            for linea in d.split("\n")
            {
                let punti: Vec<&str> = linea.split("|").collect();
                let mut punticalcolati: Vec<Point3<f64>> = Vec::new();
                
                for punto in punti
                {
                    let coords: Vec<f64> = punto.split(" ").map(|e|{e.parse::<f64>().unwrap()}).collect();
                    punticalcolati.push(
                        Point3::new(coords[0], coords[1], coords[2])
                    );
                }
                ret[i].push(
                    (
                        punticalcolati[0],
                        punticalcolati[1],
                        punticalcolati[2]
                    )
                );
            }
        }
        return ret;
    }

    /// Dalla transform, restituisce i vertici e indici
    fn plane_face_from_transform(tr: &TriLinearInterpolation, starting_index: i32) -> (Vec<Point3<f64>>,  Vec<i32>)
    {
        return (
            vec![
                tr.transform(0.0, 0.0, 0.0),
                tr.transform(1.0, 0.0, 0.0),
                tr.transform(1.0, 1.0, 0.0),
                tr.transform(0.0, 1.0, 0.0)
            ],
            vec![
                starting_index + 0,
                starting_index + 3,
                starting_index + 2,
                starting_index + 0,
                starting_index + 2,
                starting_index + 1
            ]
        );
    }

    pub trait StitchFace
    {
        /// Vertici fatti creando nuovi punti a ogni faccia
        /// Non bellissimo ma comodo per le UV
        fn vertex_coords(&self, starting_index: i32) -> Vec<Point3<f64>>
        {
            self.all_vertex_data(starting_index).0
        }
        fn vertex_indices(&self, starting_index: i32) -> Vec<i32>
        {
            self.all_vertex_data(starting_index).1
        }
        fn all_vertex_data(&self, starting_index: i32) -> (Vec<Point3<f64>>,  Vec<i32>);

        /// Serve per una visualizzazione della struttura mesh calcolata
        /// 
        /// Restituisce ogni faccia come quadrata o triangolare, invece di generare
        /// le maglie tridimensionali vere e proprie
        fn mesh_only_vertex_data(&self, starting_index: i32) -> (Vec<Point3<f64>>,  Vec<i32>);

        /// Restituisce tutti i vertici di confine della faccia,
        /// senza triangolarli.
        /// 
        /// Sono restituiti in senso antiorario
        fn boundary_vertices(&self) -> [Vec<Point3<f64>>;2];


        /// Mi dice i vertici di ogni faccia, utili per interpolare
        /// 
        /// Restituisce il tipo, i vertici interni e i vertici esterni
        fn face_data(&self) -> Vec<(StitchType, Vec<Point3<f64>>, Vec<Point3<f64>>)>;

        fn generate_transformed_vertices(&self,
            tr: &TriLinearInterpolation,
            punti_curve: &Vec<Vec<(Point3<f64>, Point3<f64>,Point3<f64>)>>,
            starting_index: i32
        ) -> (Vec<Point3<f64>>,  Vec<i32>)
        {
            let mut curve: Vec<Cylinder> = Vec::new();
            for singlecurvedata in punti_curve.iter()
            {
                for i in 1..singlecurvedata.len()
                {
                    let partenza = singlecurvedata[i-1].1;
                    let arrivo = singlecurvedata[i].1;
                    let der1 = singlecurvedata[i-1].2;
                    let der2 = singlecurvedata[i].0;
                    let bez = Bezier::new(
                        tr.transform_point(partenza),
                        tr.transform_point(arrivo),
                        Some(
                            match i==1
                            {
                                true => {tr.correct_control_point(der1, partenza)},
                                false => {tr.transform_point(der1)}
                            }
                            
                        ),
                        Some(
                            match i==singlecurvedata.len()-1
                            {
                                true => {tr.correct_control_point(der2, arrivo)},
                                false => {tr.transform_point(der2)}
                            }
                        )
                    );
                    curve.push(
                        Cylinder::new_no_taulut(Box::new(bez), get_yarn_radius(),get_yarn_radius())
                    );
                }
            }

            let mut vertici: Vec<Point3<f64>> = Vec::new();
            let mut indici: Vec<i32> = Vec::new();

            let mut last_start_index = starting_index;
            for c in curve
            {
                let res = c.to_triangles(
                    last_start_index
                );
                last_start_index += res.0.iter().max().unwrap() - last_start_index + 1;
                indici.extend(res.0);
                vertici.extend(res.1);
            }
            
            return (vertici, indici);
        }



        fn get_stitch_type(&self) -> StitchType;
    }

    impl std::fmt::Debug for dyn StitchFace
    {
        fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            write!(f, "{}", "ddd")
        }
    }

    #[derive(Debug)]
    pub struct SingleCrochetFace {
        stitch_type: StitchType,
        flip: bool,
        transform: TriLinearInterpolation
    }
    impl StitchFace for SingleCrochetFace
    {
        fn get_stitch_type(&self) -> StitchType
        {
            return self.stitch_type.clone();
        }
        fn all_vertex_data(&self, starting_index: i32) -> (Vec<Point3<f64>>,  Vec<i32>)
        {
            return self.generate_transformed_vertices(&self.transform, &SINGLE_STITCH_CURVES, starting_index);
        }

        fn mesh_only_vertex_data(&self, starting_index: i32) -> (Vec<Point3<f64>>,  Vec<i32>)
        {
            let p = plane_face_from_transform(&self.transform, starting_index);

            return p;
        }

        fn boundary_vertices(&self) -> [Vec<Point3<f64>>;2]
        {
            return [
                vec![
                    self.transform.transform(0.0, 1.0, 0.5),
                    self.transform.transform(1.0, 1.0, 0.5),
                    self.transform.transform(1.0, 0.0, 0.5),
                    self.transform.transform(0.0, 0.0, 0.5),
                ],
                Vec::new()
            ];
        }

        fn face_data(&self) -> Vec<(StitchType, Vec<Point3<f64>>, Vec<Point3<f64>>)>
        {
            if self.stitch_type == StitchType::Increase
            {
                panic!("Valore non consenito");
            }
            let faces = self.transform.get_faces();
            return vec![(
                self.stitch_type,
                faces.0,
                faces.1
            )];
        }
    }
    impl SingleCrochetFace {
        pub fn new(points: [TauCoords; 4], surf: &Box<dyn UVSurface>) -> Box<dyn StitchFace> {
            return Box::new(SingleCrochetFace {
                stitch_type: StitchType::SingleCrochet,
                flip: false,
                transform: TriLinearInterpolation::new_on_surface(points, surf)
            });
        }
    }

    
    #[derive(Debug)]
    pub struct IncreaseFace {
        stitch_type: StitchType,
        flip: bool,
        transform: [TriLinearInterpolation; 2]
    }
    impl StitchFace for IncreaseFace
    {
        fn get_stitch_type(&self) -> StitchType
        {
            return self.stitch_type.clone();
        }
        fn all_vertex_data(&self, starting_index: i32) -> (Vec<Point3<f64>>,  Vec<i32>)
        {
            let mut p1 = self.generate_transformed_vertices(&self.transform[0], &INCREASE_L_STITCH_CURVES, starting_index);
            let p2 = self.generate_transformed_vertices(&self.transform[1], &INCREASE_R_STITCH_CURVES, starting_index + p1.0.len() as i32);
            p1.0.extend(p2.0);
            p1.1.extend(p2.1);
            return p1;
        }

        fn mesh_only_vertex_data(&self, starting_index: i32) -> (Vec<Point3<f64>>,  Vec<i32>)
        {
            let mut p1 = plane_face_from_transform(&self.transform[0], starting_index);
            let mut p2 = plane_face_from_transform(&self.transform[1], starting_index + p1.0.len() as i32);

            p1.0.extend(p2.0);
            p1.1.extend(p2.1);

            return p1;
        }

        fn boundary_vertices(&self) -> [Vec<Point3<f64>>; 2]
        {
            if crochet_consts::crochet_consts::get_triangulate_n_gon()
            {
                return [
                    vec![
                        self.transform[0].transform(0.0, 0.0, 0.5),
                        self.transform[0].transform(0.0, 1.0, 0.5),
                        self.transform[0].transform(1.0, 1.0, 0.5),
                    ],
                    vec![
                        self.transform[0].transform(0.0, 0.0, 0.5),
                        self.transform[1].transform(0.0, 1.0, 0.5),
                        self.transform[1].transform(1.0, 1.0, 0.5),
                        self.transform[1].transform(1.0, 0.0, 0.5),
                        
                    ]
                ];
            }
            else
            {
                return [
                    vec![
                        self.transform[0].transform(0.0, 0.0, 0.5),
                        self.transform[0].transform(0.0, 1.0, 0.5),
                        self.transform[1].transform(0.0, 1.0, 0.5),
                        self.transform[1].transform(1.0, 1.0, 0.5),
                        self.transform[1].transform(1.0, 0.0, 0.5),
                       
                    ],
                    Vec::new()
                ];
            }
            
        }

        fn face_data(&self) -> Vec<(StitchType, Vec<Point3<f64>>, Vec<Point3<f64>>)>
        {
            #[cfg(debug_assertions)]
            {
                println!("Increase");
            }
            
            let faces_prev = self.transform[0].get_faces();
            let faces_next = self.transform[1].get_faces();
            let ret = vec![
                (
                    self.stitch_type,
                    vec![
                        faces_prev.0[0],
                        faces_next.0[1],
                        faces_next.0[2],
                        faces_next.0[3],
                        faces_prev.0[3],
                    ],
                    vec![
                        faces_prev.1[0],
                        faces_next.1[1],
                        faces_next.1[2],
                        faces_next.1[3],
                        faces_prev.1[3],
                    ]
                )
            ];
            //println!("{:?}", ret);
            return ret;
        }
    }
    impl IncreaseFace {
        pub fn new(points: [TauCoords; 5], surf: &Box<dyn UVSurface>) -> Box<dyn StitchFace> {
            let middle = points[0].lerp(&points[1], 0.5);
            return Box::new(IncreaseFace {
                stitch_type: StitchType::Increase,
                flip: false,
                transform: [
                    TriLinearInterpolation::new_on_surface([points[0], middle, points[3],points[4]], surf),
                    TriLinearInterpolation::new_on_surface([middle, points[1], points[2], points[3]], surf)
                    ]
            });
        }
    }

    #[derive(Debug)]
    pub struct DecreaseFace {
        stitch_type: StitchType,
        flip: bool,
        transform: [TriLinearInterpolation; 2]
    }
    impl StitchFace for DecreaseFace
    {
        fn get_stitch_type(&self) -> StitchType
        {
            return self.stitch_type.clone();
        }
        fn all_vertex_data(&self, starting_index: i32) -> (Vec<Point3<f64>>,  Vec<i32>)
        {
            let mut p1 = self.generate_transformed_vertices(&self.transform[0], &DECREASE_L_STITCH_CURVES, starting_index);
            let p2 = self.generate_transformed_vertices(&self.transform[1], &DECREASE_R_STITCH_CURVES, starting_index + p1.0.len() as i32);
            p1.0.extend(p2.0);
            p1.1.extend(p2.1);
            return p1;
        }

        fn mesh_only_vertex_data(&self, starting_index: i32) -> (Vec<Point3<f64>>,  Vec<i32>)
        {
            let mut p1 = plane_face_from_transform(&self.transform[0], starting_index);
            let mut p2 = plane_face_from_transform(&self.transform[1], starting_index + p1.0.len() as i32);

            p1.0.extend(p2.0);
            p1.1.extend(p2.1);

            return p1;
        }

        fn boundary_vertices(&self) -> [Vec<Point3<f64>>; 2]
        {
            if crochet_consts::crochet_consts::get_triangulate_n_gon()
            {
                return [
                    vec![
                        self.transform[0].transform(0.0, 0.0, 0.5),
                        self.transform[0].transform(0.0, 1.0, 0.5),
                        self.transform[0].transform(1.0, 0.0, 0.5),
                    ],
                    vec![
                        self.transform[1].transform(0.0, 0.0, 0.5),
                        self.transform[0].transform(0.0, 1.0, 0.5),
                        self.transform[1].transform(1.0, 1.0, 0.5),
                        self.transform[1].transform(1.0, 0.0, 0.5),
                    ]
                ];
            }
            else
            {
                return [
                    vec![
                        self.transform[0].transform(0.0, 0.0, 0.5),
                        self.transform[0].transform(0.0, 1.0, 0.5),
                        self.transform[1].transform(1.0, 1.0, 0.5),
                        self.transform[1].transform(1.0, 0.0, 0.5),
                        self.transform[0].transform(1.0, 0.0, 0.5),
                    ],
                    Vec::new()
                ];
            }
            
        }

        fn face_data(&self) -> Vec<(StitchType, Vec<Point3<f64>>, Vec<Point3<f64>>)>
        {
            #[cfg(debug_assertions)]
            {
                println!("Decrease");
            }
            let faces_prev = self.transform[0].get_faces();
            let faces_next = self.transform[1].get_faces();
            let ret = vec![
                (
                    self.stitch_type,
                    vec![
                        faces_prev.0[0],
                        faces_next.0[1],
                        faces_next.0[2],
                        faces_next.0[3],
                        faces_prev.0[3],
                    ],
                    vec![
                        faces_prev.1[0],
                        faces_next.1[1],
                        faces_next.1[2],
                        faces_next.1[3],
                        faces_prev.1[3],
                    ]
                )
            ];
            return ret;
        }
    }
    impl DecreaseFace {
        pub fn new(points: [TauCoords; 5], surf: &Box<dyn UVSurface>) -> Box<dyn StitchFace> {
            let middle = points[3].lerp(&points[4], 0.5);
            return Box::new(DecreaseFace {
                stitch_type: StitchType::Decrease,
                flip: false,
                transform: [
                    TriLinearInterpolation::new_on_surface([points[0], points[1], middle, points[4]], surf),
                    TriLinearInterpolation::new_on_surface([points[1], points[2], points[3], middle], surf)
                    ]
            });
        }
    }
    

    #[derive(Debug)]
    pub struct TurnBackFace {
        stitch_type: StitchType,
        flip: bool,
        transform: TriLinearInterpolation
    }
    impl StitchFace for TurnBackFace
    {
        fn get_stitch_type(&self) -> StitchType
        {
            return self.stitch_type.clone();
        }
        fn all_vertex_data(&self, starting_index: i32) -> (Vec<Point3<f64>>,  Vec<i32>)
        {
            let p1 = self.generate_transformed_vertices(&self.transform, &TURNBACK_STITCH_CURVES, starting_index);
            return p1;
        }

        fn mesh_only_vertex_data(&self, starting_index: i32) -> (Vec<Point3<f64>>,  Vec<i32>)
        {
            let p = plane_face_from_transform(&self.transform, starting_index);

            return p;
        }

        fn boundary_vertices(&self) -> [Vec<Point3<f64>>; 2]
        {
            return [
                vec![
                    self.transform.transform(0.0, 0.0, 0.5),
                    self.transform.transform(1.0, 0.0, 0.5),
                    self.transform.transform(1.0, 1.0, 0.5),
                    self.transform.transform(0.0, 1.0, 0.5)
                ],
                Vec::new()
            ];
        }

        fn face_data(&self) -> Vec<(StitchType, Vec<Point3<f64>>, Vec<Point3<f64>>)>
        {
            let faces = self.transform.get_faces();
            return vec![(
                self.stitch_type,
                faces.0,
                faces.1
            )];
        }
    }
    impl TurnBackFace {
        pub fn new(points: [TauCoords; 4], surf: &Box<dyn UVSurface>) -> Box<dyn StitchFace> {
            return Box::new(TurnBackFace {
                stitch_type: StitchType::TurnBack,
                flip: false,
                transform: 
                    TriLinearInterpolation::new_on_surface([points[0], points[1], points[2], points[3]], surf)
            });
        }
    }

    #[derive(Debug)]
    pub struct TurnForwardFace {
        stitch_type: StitchType,
        flip: bool,
        transform: TriLinearInterpolation
    }
    impl StitchFace for TurnForwardFace
    {
        fn get_stitch_type(&self) -> StitchType
        {
            return self.stitch_type.clone();
        }
        fn all_vertex_data(&self, starting_index: i32) -> (Vec<Point3<f64>>,  Vec<i32>)
        {
            let p1 = self.generate_transformed_vertices(&self.transform, &TURNFORWARD_STITCH_CURVES, starting_index);
            return p1;
        }

        fn mesh_only_vertex_data(&self, starting_index: i32) -> (Vec<Point3<f64>>,  Vec<i32>)
        {
            let p = plane_face_from_transform(&self.transform, starting_index);

            return p;
        }

        fn boundary_vertices(&self) -> [Vec<Point3<f64>>; 2]
        {
            return [
                vec![
                    self.transform.transform(0.0, 0.0, 0.5),
                    self.transform.transform(1.0, 0.0, 0.5),
                    self.transform.transform(1.0, 1.0, 0.5),
                    self.transform.transform(0.0, 1.0, 0.5)
                ],
                Vec::new()
            ];
        }

        fn face_data(&self) -> Vec<(StitchType, Vec<Point3<f64>>, Vec<Point3<f64>>)>
        {
            let faces = self.transform.get_faces();
            return vec![(
                self.stitch_type,
                faces.0,
                faces.1
            )];
        }
    }
    impl TurnForwardFace {
        pub fn new(points: [TauCoords; 4], surf: &Box<dyn UVSurface>) -> Box<dyn StitchFace> {
            return Box::new(TurnForwardFace {
                stitch_type: StitchType::TurnForward,
                flip: false,
                transform: 
                    TriLinearInterpolation::new_on_surface([points[0], points[1], points[2], points[3]], surf)
            });
        }
    }




    #[derive(Debug)]
    pub struct CapFace {
        stitch_type: StitchType,
        flip: bool,
        transform: TriLinearInterpolation
    }
    impl StitchFace for CapFace
    {
        fn get_stitch_type(&self) -> StitchType
        {
            return self.stitch_type.clone();
        }
        fn all_vertex_data(&self, starting_index: i32) -> (Vec<Point3<f64>>,  Vec<i32>)
        {
            return self.generate_transformed_vertices(&self.transform, &CAP_STITCH_CURVES, starting_index);
        }

        fn mesh_only_vertex_data(&self, starting_index: i32) -> (Vec<Point3<f64>>,  Vec<i32>)
        {
            let p = plane_face_from_transform(&self.transform, starting_index);

            return p;
        }

        fn boundary_vertices(&self) -> [Vec<Point3<f64>>; 2]
        {
            return [Vec::new(), Vec::new()];
        }

        fn face_data(&self) -> Vec<(StitchType, Vec<Point3<f64>>, Vec<Point3<f64>>)>
        {
            let faces = self.transform.get_faces();
            return vec![(
                self.stitch_type,
                faces.0,
                faces.1
            )];
        }
    }
    impl CapFace {
        pub fn new(points_internal: [Point3<f64>; 4], points_ext: [Point3<f64>; 4]) -> Box<dyn StitchFace> {
            return Box::new(CapFace {
                stitch_type: StitchType::EndingCap,
                flip: false,
                transform: TriLinearInterpolation::new(points_internal, points_ext)
            });
        }
    }
}
