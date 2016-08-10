use vectorspace::{VectorSpace};
use bspline::{Bspline,SplineData,KnotManip,ClassInvariant, Curve};

pub struct RationalBspline<Point:VectorSpace>  {
    spl: Bspline<Point::H>,
}

impl <Point> RationalBspline<Point> where Point:VectorSpace
{
    pub fn new(control_points: Vec<Point::H>,  knots: Vec<f64>) -> RationalBspline<Point>
    {
        let spl = Bspline::new( control_points, knots);
        RationalBspline{spl : spl}
    }
}

impl <Point> SplineData for RationalBspline<Point> where Point:VectorSpace
{
    type T = Point::H;
    fn control_points(&self) -> &Vec<Point::H> { self.spl.control_points() }
    fn knots(&self) -> &Vec<f64>  { self.spl.knots() }
    fn degree(&self) -> u32 { self.spl.degree() }
}

impl<P> ClassInvariant for RationalBspline<P> where P:VectorSpace
{
    fn is_valid(&self) -> Result<bool, &str> {
        self.spl.is_valid()
    }
}

// inputs are of form v,v',v"..v^(n) (ie all derivatives upto n-th)
// returns w^(n)
fn rational_derivatives_from_derivatives<T>(vecs: &Vec<T>)-> Vec<T::L>
    where T:VectorSpace
{
    let mut vs : Vec<T::L> = Vec::with_capacity(vecs.len());
    let der_order = vecs.len();
    let mut bbasis:Vec<f64> = vec![0.0;der_order];
    let dim = T::dim() - 1;
    let mut ws:Vec<f64> = vec![1.0;der_order];
    for i in 0..der_order
    {
        for j in (1..i).rev() {
            bbasis[j]+=bbasis[j-1];
        }
        bbasis[i] = 1.0;
        let mut v = vecs[i].ldim();
        let w = vecs[i].extract(dim);
        ws.push(w);
        for j in 1..i+1 {
            let mut u = vs[i-j];
            u = u*bbasis[j]*ws[j];
            v = v - u
        }
        v = v*(1.0/vecs[0].extract(dim));
        vs.push(v);
    }
    vs
}

impl<P> Curve for RationalBspline<P> where P:VectorSpace
{
    type T = <<P as VectorSpace>::H as VectorSpace>::L;
    fn param_range(&self)->(f64,f64) { self.spl.param_range() }
    fn eval(&self,v:f64) -> Self::T {
        let pw = self.spl.eval(v);
        pw.ldim()*(1.0/pw.extract(P::dim()))
    }
    fn eval_derivative(&self, v:f64,order:u32) -> Self::T
    {
        let mut vecs : Vec<P::H> = Vec::with_capacity(order as usize+1);
        for i in 0..order+1 {
            vecs.push(self.spl.eval_derivative(v,i))
        }
        rational_derivatives_from_derivatives(&vecs)[0]
    }
}

impl<T> KnotManip for RationalBspline<T> where T:VectorSpace {
    fn start_mult(&self) -> usize {
        self.spl.start_mult()
    }

    fn end_mult(&self) -> usize {
        self.spl.end_mult()
    }

    fn front(&self) -> f64 {
        self.spl.front()
    }

    fn back(&self) -> f64 {
        self.spl.back()
    }

    fn locate_nu(&self, u: f64) -> usize {
        self.spl.locate_nu(u)
    }

    fn mult(&self, u: f64) -> usize {
        self.spl.mult(u)
    }

    fn rebase(&self, taus: Vec<f64>) -> RationalBspline<T> {
        RationalBspline{ spl : self.spl.rebase(taus) }
    }

    fn insert_knot(&self, tau: f64) -> RationalBspline<T> {
        RationalBspline{ spl : self.spl.insert_knot(tau) }
    }

    fn insert_knots(&self, taus: &Vec<f64>) -> RationalBspline<T> {
        RationalBspline{ spl : self.spl.insert_knots(taus) }
    }
}
#[test]
fn it_works() {
    use point::{Point2,Point1};
    let v : Vec<Point2> = vec![Point2::new(1.0,1.0), Point2::new(2.0,1.0),Point2::new(1.0, 1.0)];
    type RsP1 = RationalBspline<Point1>;
    let rs = RsP1::new(v, vec!(0.,0.,0.,1.,1.,1.) );

    let bs = Bspline::new(vec![1.0,2.0,1.0], vec![1.,1.,1.,2.,2.,2.]);
    assert_eq!(rs.eval(0.0) , bs.eval(1.0));
    // assert_eq!((rs.eval_derivative(0.0,1)-bs.eval_derivative(1.0,1)).len(), 0.0);
}
