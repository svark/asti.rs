use vectorspace::VectorSpace;
pub trait Curve
{
    type T:VectorSpace;
    fn param_range(&self)->(f64,f64);
    fn eval(&self,v:f64) -> Self::T;
    fn eval_derivative(&self, v:f64,order:u32) -> Self::T;
}

pub trait BlossomCurve
{
    type T : VectorSpace;
    fn blossom_eval(&self, der_order: u32, us: &[f64]) -> Self::T;
}
