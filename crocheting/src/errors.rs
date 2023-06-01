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

pub mod errors
{
    use std::fmt;

    #[derive(Debug)]
    #[derive(PartialEq)]
    #[derive(Clone)]
    pub enum ErrorType
    {
        ZeroVector,
        RadiusInfiniteDistance,
        OutOfRange,// Per esempio se chiedo un "u" fuori da 0-1
    }

    #[derive(Debug)]
    #[derive(Clone)]
    pub struct CrochetError
    {
        pub description: String,
        pub error_type: ErrorType
    }

    impl std::fmt::Display for CrochetError {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "{}", self.description)
        }
    }
}