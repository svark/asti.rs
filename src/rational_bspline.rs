use vectorspace::VectorSpace;
use bspline::{Bspline, SplineWrapper, SplineMut};
use periodic_spline::PeriodicBspline;
use curve::Curve;
use class_invariant::ClassInvariant;

pub struct RationalBspline<Point>
    where Point: VectorSpace
{
    spl: Bspline<Point::H>,
}

pub struct PeriodicRationalBspline<Point>
    where Point: VectorSpace
{
    spl: PeriodicBspline<Point::H>,
}


pub trait RationalTrait: SplineWrapper
{
    // inputs are of form v,v',v"..v^(n) (ie all derivatives upto n-th)
    // returns w^(n)
    fn rational_derivatives_from_derivatives(vecs: &Vec< <Self as SplineWrapper>::TW>) -> Vec<<<Self as SplineWrapper>::TW as VectorSpace>::L>
    {
        let mut vs: Vec<<<Self as SplineWrapper>::TW as VectorSpace>::L> = Vec::with_capacity(vecs.len());
        let der_order = vecs.len();
        let mut bbasis: Vec<f64> = vec![0.0;der_order];
        let dim = Self::TW::dim() - 1;
        let mut ws: Vec<f64> = Vec::with_capacity(der_order + 1);
        for i in 0..der_order {
            for j in (1..i).rev() {
                bbasis[j] += bbasis[j - 1];
            }
            bbasis[i] = 1.0;
            let mut v = vecs[i].ldim();
            let w = vecs[i].extract(dim);
            ws.push(w);
            for j in 1..i + 1 {
                let mut u = vs[i - j];
                u = u * bbasis[j] * ws[j];
                v = v - u
            }
            v = v * (1.0 / vecs[0].extract(dim));
            vs.push(v);
        }
        vs
    }

}


impl<Point> SplineWrapper for RationalBspline<Point>
    where Point: VectorSpace
{
    type TW = Point::H;
    fn to_spline(&self) -> &Bspline<Self::TW> {
        &self.spl
    }
    fn from_spline(spl: Bspline<Self::TW>) -> Self {
        RationalBspline { spl: spl }
    }
}

impl<Point> SplineWrapper for PeriodicRationalBspline<Point>
    where Point: VectorSpace
{
    type TW = Point::H;
    fn to_spline(&self) -> &Bspline<Self::TW> {
        self.spl.to_spline()
    }
    fn from_spline(spl: Bspline<Self::TW>) -> Self {
        PeriodicRationalBspline { spl: PeriodicBspline::from_spline(spl) }
    }
}

impl<Point: VectorSpace> RationalTrait for RationalBspline<Point> {}

impl<Point: VectorSpace> RationalTrait for PeriodicRationalBspline<Point> {}


impl<P> ClassInvariant for RationalBspline<P>
    where P: VectorSpace
{
    fn is_valid(&self) -> Result<bool, &str> {
        self.spl.is_valid()
    }
}

impl<SplineType: RationalTrait> Curve for SplineType {
    type T = <<Self as SplineWrapper>::TW as VectorSpace>::L;
    fn param_range(&self) -> (f64, f64) {
        self.to_spline().param_range()
    }

    fn eval(&self, v: f64) -> Self::T {
        let pw = self.to_spline().eval(v);
        pw.ldim() * (1.0 / pw.extract(Self::T::dim()))
    }

    fn eval_derivative(&self, v: f64, order: u32) -> Self::T {
        let mut vecs: Vec<<Self as SplineWrapper>::TW> = Vec::with_capacity(order as usize + 1);
        for i in 0..order + 1 {
            vecs.push(self.to_spline().eval_derivative(v, i))
        }
        *SplineType::rational_derivatives_from_derivatives(&vecs).last().unwrap()
    }
}

impl<P: VectorSpace> SplineMut for RationalBspline<P> {
    fn into_spline(self) -> Bspline<Self::T> {
        self.spl
    }
}

impl<P: VectorSpace> SplineMut for PeriodicRationalBspline<P> {
    fn into_spline(self) -> Bspline<Self::T> {
        self.spl.into_spline()
    }
}



#[test]
fn it_works() {
    use point::{Point2, Point1};
    use splinedata::SplineData;
    let v: Vec<Point2> = vec![Point2::new(1.0, 1.0), Point2::new(2.0, 1.0), Point2::new(1.0, 1.0)];
    type RsP1 = RationalBspline<Point1>;
    let rs = RsP1::new(v, vec![0., 0., 0., 1., 1., 1.]);

    let bs = Bspline::new(vec![1.0, 2.0, 1.0], vec![0., 0.0, 0.0, 1., 1., 1.]);
    assert_eq!(rs.eval(0.0), bs.eval(1.0));
    assert_eq!(rs.eval_derivative(0.0, 1), bs.eval_derivative(0.0, 1));
}
