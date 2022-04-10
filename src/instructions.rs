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

pub mod instructions
{
    use crate::mesh_viz::mesh_viz::{StitchMesh, StitchConnection, StitchType};
    use std::fs::File;
    use std::io::Write;


    pub struct InstructionsGenerator
    {
        computed_instructions: Vec<Vec<(StitchType, u32)>>
    }
    impl InstructionsGenerator
    {
        pub fn new(st: &StitchMesh) -> Self
        {
            // Logica creazione istruzioni
            let tracewale = st.get_trace(StitchConnection::NextWale);
            let tracecourse = st.get_trace(StitchConnection::NextCourse);
            let mut istr: Vec<Vec<(StitchType, u32)>> = Vec::new();

            {
                // Separo in vettori per ogni riga
                let mut w_ind: usize = 0;

                for cur_instr_data in tracecourse.iter()
                {
                    let cur_start_row = tracewale[w_ind % tracewale.len()];
                    if cur_start_row.0 == cur_instr_data.0
                    {
                        w_ind += 1;
                        istr.push(Vec::new());
                    }
                    istr.last_mut().unwrap().push(
                        (
                            cur_instr_data.1,
                            1
                        )
                    );
                }
            }
            

            // Collassa
            for row_index in 0..istr.len()
            {
                let mut inner_index = 1;
                while inner_index < istr[row_index].len()
                {
                    if istr[row_index][inner_index].0 == istr[row_index][inner_index-1].0
                    {
                        // Allora le collasso
                        // Rimuovo il secondo e lo metto nel primo
                        let removed = istr[row_index].remove(inner_index);
                        istr[row_index][inner_index-1].1 += removed.1;
                    }
                    else
                    {
                        inner_index += 1;
                    }
                }
            }
            
            //
            return InstructionsGenerator{computed_instructions: istr}
        }

        pub fn text_instructions(&self, filename: Option<&str>) -> String
        {
            #[derive(PartialEq)]
            enum KnittingDirection
            {
                Forward,
                Back
            }

            let mut outp = String::new();

            let mut direzione: KnittingDirection = KnittingDirection::Forward;

            let mut catenelle = 0;
            for istr in self.computed_instructions[0].iter()
            {
                let v = match istr.0
                {
                    StitchType::SingleCrochet => 1,
                    StitchType::Increase => 1,
                    StitchType::Decrease => 2,
                    StitchType::SlipSitch => 1,
                    _ => 0
                };
                catenelle += v*istr.1;
            }
            outp.push_str(&format!(
                "Anello magico:\n    {} catenelle\n\n", catenelle
            ));

            for (row_index, row_data) in self.computed_instructions.iter().enumerate()
            {
                outp.push_str(&format!(
                    "Giro {}\n", row_index+1
                ));

                for stitch in row_data.iter()
                {
                    outp.push_str(&format!(
                            "    {}{} {}\n",
                            match direzione
                            {
                                KnittingDirection::Back => {
                                    if stitch.0 != StitchType::TurnBack
                                    &&
                                    stitch.0 != StitchType::TurnForward
                                    {
                                        "-"
                                    }
                                    else
                                    {
                                        ""
                                    }
                                },
                                _ => ""
                            },
                            // numero
                            match stitch.0 {
                                StitchType::TurnBack | StitchType::TurnForward => String::new(),
                                _ => stitch.1.to_string()
                            },
                            match stitch.0 {
                                StitchType::SingleCrochet => "mb",
                                StitchType::Decrease => "dim",
                                StitchType::Increase => "aum",
                                StitchType::TurnBack => {
                                    direzione = KnittingDirection::Back;
                                    "←"
                                },
                                StitchType::TurnForward => {
                                    direzione = KnittingDirection::Forward;
                                    "→"
                                }
                                _ => panic!("Non supportato")
                            }
                        )
                    );
                }
                outp.push_str("\n");
            }
            
            if filename.is_some()
            {
                let mut output = File::create(filename.unwrap()).unwrap();
                write!(output, "{}", outp).unwrap();
            }
            return outp;
        }
    }
}