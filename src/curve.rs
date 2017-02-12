use vectorspace::PointT;
pub trait Curve
{
    type T:PointT;
    fn eval(&self, v: f64) -> Self::T;
    fn eval_derivative(&self, v: f64, order: u32) -> Self::T;
}

pub trait Domain
{
    fn param_range(&self) -> (f64, f64);
    fn start_param(&self) -> f64 {
        self.param_range().0
    }
    fn end_param(&self) -> f64 {
        self.param_range().1
    }
}

pub trait FiniteCurve: Curve + Domain
{

}

impl<T> FiniteCurve for T where T: Domain + Curve {}

pub trait BlossomCurve
{
    type T : PointT;
    fn blossom_eval(&self, der_order: u32, us: &[f64]) -> Self::T;
}
// todo use lazy thunk
#[derive(Debug,Copy,Clone)]
pub struct CurvePoint<P: PointT> {
    pub pnt: P,
    pub t: f64,
}
