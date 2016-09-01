use vectorspace::VectorSpace;

pub trait SplineData
{
    type T : VectorSpace;
    fn control_points(&self) -> &Vec<Self::T>;
    fn knots(&self) -> &Vec<f64> ;
    fn degree(&self) -> u32;
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


pub fn  greville<SplineType>(spl : &SplineType, k : usize ) -> f64
    where SplineType:SplineData
{
    let d = spl.degree() as usize;
    let mut  sum = 0.0;
    sum = spl.knots()[k+1..k+d+1].iter().fold(sum, |s,x| { s + x });
    sum/(d as f64)
}
