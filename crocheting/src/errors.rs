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