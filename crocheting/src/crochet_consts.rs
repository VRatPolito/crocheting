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

pub mod crochet_consts
{
    pub static mut STITCH_WIDTH: f64 = 0.01;
    pub static mut STITCH_HEIGHT: f64 = 0.01;
    pub static mut GENERATE_SHORT_ROWS: bool = true;

    pub const ZERO_VECTOR_TOLERANCE: f64 = 0.00001;
    pub const ZERO_CURVATURE_TOLERANCE: f64 = 0.00001;

    pub static mut SNAP_TRESHOLD: f64 = 0.3;

    // Le facce di aumento e diminuzione sono facce da cinque vertici
    // Per coerenza con il paper di Capunaman, possono essere divise
    // in un triangolo e un quadrilatero
    pub static mut TRIANGULATE_N_GON: bool = true;


    // Usare meccanismo a doppia soglia?
    pub static mut USE_DOUBLE_THRESHOLD: bool = true;


    // Da che valore di V incomincia la spirale?
    pub static mut V_ROTATION_SURFACE: f64 = 0.0;


    pub static mut DIAGONAL_HEIGHT_INCREMENT: bool = false;


    // Variabili per sensibilità divisione e giunzione maglie
    const SPLIT_THRESH: f64 = 1.25;
    const JOIN_THRESH: f64 = 0.75;
    const SPLIT_THRESH_SHRINK: f64 = 1.5;
    const JOIN_THRESH_SHRINK: f64 = 1.0;
    const SPLIT_THRESH_EXPANDING: f64 = 1.0;
    const JOIN_THRESH_EXPANDING: f64 = 0.5;
    

    // Fattore amplificazione segnale dithering
    const DITHERING_MULT: f64 = 0.15;
    const DITHERING_MAX_TICK: usize = 8;
    static mut DITHERING_CUR_TICK: usize = 0;

    // I due numeri qua sotto non sono identici
    // per creare una isteresi e renderlo più stabile
    /// Cross product deve arrivare a -1 per la contrazione, +1 per l'espansione
    pub fn target_to_split_width(cross_product: f64) -> f64
    {
        // return SPLIT_THRESH;
        // return 2.5;
        let p = cross_product.clamp(-1.0, 1.0) / 4.0;
        if p > 0.0
        {
            return SPLIT_THRESH*(1.0-p) + SPLIT_THRESH_EXPANDING*p;
        }
        else
        {
            let p = p.abs();
            return SPLIT_THRESH*(1.0-p) + SPLIT_THRESH_SHRINK*p;
        }
    }
    pub fn target_to_join_width(cross_product: f64) -> f64
    {
        // return JOIN_THRESH;
        // return 2.5;
        let p = cross_product.clamp(-1.0, 1.0) / 4.0;
        if p > 0.0
        {
            return JOIN_THRESH*(1.0-p) + JOIN_THRESH_EXPANDING*p;
        }
        else
        {
            let p = p.abs();
            return JOIN_THRESH*(1.0-p) + JOIN_THRESH_SHRINK*p;
        }
    }

    pub fn get_width() -> f64
    {
        unsafe
        {
            return STITCH_WIDTH.clone()
        }
    }
    pub fn get_height() -> f64
    {
        unsafe
        {
            return STITCH_HEIGHT.clone()
        }
    }

    pub fn get_thickness() -> f64
    {
        unsafe
        {
            return STITCH_HEIGHT / 2.0;
        }
    }

    pub fn get_yarn_radius() -> f64
    {
        match get_width() < get_height()
        {
            true => return get_width()*0.05,
            false => return get_height()*0.05
        }
    }

    /// Onda spaziale con frequenza 6
    /// 
    /// Ampiezza nell'intervallo 0-1
    pub fn dithering(v: f64) -> f64
    {
        // let v = unsafe{DITHERING_CUR_TICK} as f64 / DITHERING_MAX_TICK as f64;
        // return 0.0;
        const FREQUENCY: f64 = 6.0;
        let sin_value = (v * 2.0* std::f64::consts::PI * FREQUENCY)
        .sin();
        return sin_value*DITHERING_MULT;
        if sin_value > 0.0
        {
            return sin_value*DITHERING_MULT;
        }
        else
        {
            return sin_value*(DITHERING_MULT*0.5);
        }
    }
    pub fn dithering_advance_tick()
    {
        unsafe {
            DITHERING_CUR_TICK = (DITHERING_CUR_TICK + 1) % DITHERING_MAX_TICK;
        }
    }
    pub fn dithering_reset_tick()
    {
        unsafe {
            DITHERING_CUR_TICK = 0;
        }
    }

    /// Cambia la dimensione in *metri* della maglia bassa
    pub unsafe fn set_stitch_dimensions(new_width: f64, new_height: f64)
    {
        STITCH_WIDTH = new_width;
        STITCH_HEIGHT = new_height;
    }

    pub unsafe fn set_triangulate_n_gon(v: bool)
    {
        TRIANGULATE_N_GON = v;
    }
    pub fn get_triangulate_n_gon() -> bool
    {
        unsafe {return TRIANGULATE_N_GON}
    }

    pub unsafe fn set_double_threshold(v: bool)
    {
        USE_DOUBLE_THRESHOLD = v;
    }
    pub fn get_double_threshold() -> bool
    {
        unsafe {return USE_DOUBLE_THRESHOLD}
    }


    pub unsafe fn set_flag_generate_short_row_flag(new_val: bool)
    {
        GENERATE_SHORT_ROWS = new_val;
    }
    pub fn generate_short_row_flag() -> bool
    {
        unsafe
        {
            return GENERATE_SHORT_ROWS.clone();
        }
    }

    pub unsafe fn set_snap_treshold(new_val: f64)
    {
        assert!(new_val>=0.0&&new_val<=1.0);
        SNAP_TRESHOLD = new_val;
    }
    pub fn get_snap_treshold() -> f64
    {
        unsafe
        {
            return SNAP_TRESHOLD.clone();
        }
    }
    pub unsafe fn reset_snap_treshold()
    {
        SNAP_TRESHOLD = 0.0;
    }


    pub fn get_v_rotation_surface() -> f64
    {
        unsafe
        {
            return V_ROTATION_SURFACE.clone();
        }
    }
    pub unsafe fn set_v_rotation_surface(new_val: f64)
    {
        V_ROTATION_SURFACE = new_val;
    }


    pub fn get_diagonal_height_increment() -> bool
    {
        unsafe
        {
            return DIAGONAL_HEIGHT_INCREMENT.clone();
        }
    }
    pub unsafe fn set_diagonal_height_increment(new_val: bool)
    {
        DIAGONAL_HEIGHT_INCREMENT = new_val;
    }
}