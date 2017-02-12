use vectorspace::PointT;
use bspline::{Bspline, SplineWrapper};
use periodic_spline::PeriodicBspline;
use curve::Curve;
use class_invariant::ClassInvariant;
use vectorspace::{Ops, AssocPoint};

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

impl<P: PointT> From<Bspline<P::H>> for RationalBspline<P> {
    fn from(spl: Bspline<P::H>) -> RationalBspline<P> {
        RationalBspline { spl: spl }
    }
}

impl<Point: PointT> AssocPoint for RationalBspline<Point> {
    type TW = Point::H;
}

impl<Point: PointT> AssocPoint for PeriodicRationalBspline<Point> {
    type TW = Point::H;
}


impl<Point> SplineWrapper for RationalBspline<Point> where Point: PointT {}

impl<Point> SplineWrapper for PeriodicRationalBspline<Point> where Point: PointT {}

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
        let pw = self.as_ref().eval(v);
        pw.ldim() * (1.0 / pw[Self::T::dim()])
    }

    fn eval_derivative(&self, v: f64, order: u32) -> Self::T {
        let vecs: Vec<_> = (0..order + 1)
                               .into_iter()
                               .map(|i| self.as_ref().eval_derivative(v, i))
                               .collect();

        *rational_derivatives_from_derivatives(&vecs).last().unwrap()
    }
}

impl<P: PointT> AsRef<Bspline<P::H>> for RationalBspline<P> {
    fn as_ref(&self) -> &Bspline<P::H> {
        &self.spl
    }
}

impl<P: PointT> AsMut<Bspline<P::H>> for RationalBspline<P> {
    fn as_mut(&mut self) -> &mut Bspline<P::H> {
        &mut self.spl
    }
}

impl<P: PointT> AsRef<PeriodicBspline<P::H>> for PeriodicRationalBspline<P> {
    fn as_ref(&self) -> &PeriodicBspline<P::H> {
        &self.spl
    }
}

impl<P: PointT> AsMut<PeriodicBspline<P::H>> for PeriodicRationalBspline<P> {
    fn as_mut(&mut self) -> &mut PeriodicBspline<P::H> {
        &mut self.spl
    }
}

impl<P: PointT> From<PeriodicBspline<P::H>> for PeriodicRationalBspline<P> {
    fn from(spl: PeriodicBspline<P::H>) -> PeriodicRationalBspline<P> {
        PeriodicRationalBspline { spl: spl }
    }
}
// conversion to and from Bspline
impl<P: PointT> AsRef<Bspline<P::H>> for PeriodicRationalBspline<P> {
    fn as_ref(&self) -> &Bspline<P::H> {
        self.spl.as_ref()
    }
}

impl<P: PointT> AsMut<Bspline<P::H>> for PeriodicRationalBspline<P> {
    fn as_mut(&mut self) -> &mut Bspline<P::H> {
        self.spl.as_mut()
    }
}

impl<P: PointT> From<Bspline<P::H>> for PeriodicRationalBspline<P> {
    fn from(spl: Bspline<P::H>) -> PeriodicRationalBspline<P> {
        PeriodicRationalBspline { spl: PeriodicBspline::from(spl) }
    }
}

impl<P: PointT> Into<Bspline<P::H>> for PeriodicRationalBspline<P> {
    fn into(self) -> Bspline<P::H> {
        self.spl.into()
    }
}

impl<P: PointT> Into<Bspline<P::H>> for RationalBspline<P> {
    fn into(self) -> Bspline<P::H> {
        self.spl
    }
}

#[test]
fn it_works() {
    use point::{Pt2, Pt1};
    let v: Vec<Pt2> = vec![Pt2::new(1.0, 1.0), Pt2::new(2.0, 1.0), Pt2::new(1.0, 1.0)];
    let rs = RationalBspline::<Pt1>::from(Bspline::new(v, vec![0., 0., 0., 1., 1., 1.]));

    let bs = Bspline::new(vec![Pt1::new(1.0), Pt1::new(2.0), Pt1::new(1.0)],
                          vec![0., 0.0, 0.0, 1., 1., 1.]);
    assert_eq!(rs.eval(0.0), bs.eval(1.0));
    assert_eq!(rs.eval_derivative(0.0, 1), bs.eval_derivative(0.0, 1));
}
