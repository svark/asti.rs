use vectorspace::VectorSpace;

pub trait SplineData
{
    type T : VectorSpace;
    fn control_points(&self) -> &Vec<Self::T>;
    fn knots(&self) -> &Vec<f64> ;
    fn degree(&self) -> u32;
    fn new(cpts: Vec<Self::T>, ts : Vec<f64>) -> Self;
}


pub trait KnotManip
{
    fn start_mult(&self) -> usize;
    fn end_mult(&self) -> usize;
    fn mult(&self, u: f64) -> usize;
    fn front(&self) -> f64;
    fn back(&self) -> f64;
    fn locate_nu(&self, u: f64) -> usize;
    fn rebase(&self, taus: Vec<f64>) -> Self;
    fn insert_knot(&self, tau: f64) -> Self;
    fn insert_knots(&self, taus: &Vec<f64>) ->Self;
}


