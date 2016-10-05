use vectorspace::PointT;
use bspline::{Bspline, SplineWrapper, SplineMut};
use periodic_spline::PeriodicBspline;
use curve::{Curve, FiniteCurve};
use class_invariant::ClassInvariant;
use vectorspace::Ops;

pub struct RationalBspline<Point>
    where Point: PointT
{
    spl: Bspline<Point::H>,
}

pub struct PeriodicRationalBspline<Point>
    where Point: PointT
{
    spl: PeriodicBspline<Point::H>,
}


// inputs are of form v,v',v"..v^(n) (ie all derivatives upto n-th)
// returns w^(n)
pub fn rational_derivatives_from_derivatives<T: PointT>(vecs: &Vec<T>) -> Vec<T::L> {
    let mut vs: Vec<T::L> = Vec::with_capacity(vecs.len());
    let der_order_1 = vecs.len();
    let mut bbasis: Vec<f64> = vec![0.0;der_order_1];
    let dim = T::dim() - 1;
    let mut ws: Vec<f64> = Vec::with_capacity(der_order_1 + 1);
    for i in 0..der_order_1 {
        for j in (1..i).rev() {
            bbasis[j] += bbasis[j - 1];
        }
        bbasis[i] = 1.0;
        let mut v = vecs[i].ldim();
        let w = vecs[i].extract(dim);
        ws.push(w);
        for j in 1..i + 1 {
            let u = vs[i - j] * bbasis[j] * ws[j];
            v += T::L::zero_pt() - u;
        }
        v = v * (1.0 / vecs[0].extract(dim));
        vs.push(v);
    }
    vs
}

pub trait RationalTrait: SplineWrapper
{

}


impl<Point> SplineWrapper for RationalBspline<Point>
    where Point: PointT
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
    where Point: PointT
{
    type TW = Point::H;
    fn to_spline(&self) -> &Bspline<Self::TW> {
        self.spl.to_spline()
    }
    fn from_spline(spl: Bspline<Self::TW>) -> Self {
        PeriodicRationalBspline { spl: PeriodicBspline::from_spline(spl) }
    }
}

impl<Point: PointT> RationalTrait for RationalBspline<Point> {}

impl<Point: PointT> RationalTrait for PeriodicRationalBspline<Point> {}


impl<P> ClassInvariant for RationalBspline<P>
    where P: PointT
{
    fn is_valid(&self) -> Result<bool, &str> {
        self.spl.is_valid()
    }
}

impl<SplineType: RationalTrait> Curve for SplineType {
    type T = TWL!();

    fn eval(&self, v: f64) -> Self::T {
        let pw = self.to_spline().eval(v);
        pw.ldim() * (1.0 / pw[Self::T::dim()])
    }

    fn eval_derivative(&self, v: f64, order: u32) -> Self::T {
        let vecs: Vec<_> = (0..order + 1)
                               .into_iter()
                               .map(|i| self.to_spline().eval_derivative(v, i))
                               .collect();

        *rational_derivatives_from_derivatives(&vecs).last().unwrap()
    }
}

impl<SplineType: RationalTrait> FiniteCurve for SplineType {
    fn param_range(&self) -> (f64, f64) {
        self.to_spline().param_range()
    }
}

impl<P: PointT> SplineMut for RationalBspline<P> {
    fn into_spline(self) -> Bspline<Self::T> {
        self.spl
    }
}

impl<P: PointT> SplineMut for PeriodicRationalBspline<P> {
    fn into_spline(self) -> Bspline<Self::T> {
        self.spl.into_spline()
    }
}


#[test]
fn it_works() {
    use point::{Pt2, Pt1};
    let v: Vec<Pt2> = vec![Pt2::new(1.0, 1.0), Pt2::new(2.0, 1.0), Pt2::new(1.0, 1.0)];
    let rs = RationalBspline::<Pt1>::from_spline(Bspline::new(v, vec![0., 0., 0., 1., 1., 1.]));

    let bs = Bspline::new(vec![Pt1::new(1.0), Pt1::new(2.0), Pt1::new(1.0)],
                          vec![0., 0.0, 0.0, 1., 1., 1.]);
    assert_eq!(rs.eval(0.0), bs.eval(1.0));
    assert_eq!(rs.eval_derivative(0.0, 1), bs.eval_derivative(0.0, 1));
}
