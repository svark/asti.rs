use splinedata::{SplineData, KnotManip};
use bspline::{Bspline, SplineWrapper, SplineMut};
use class_invariant::ClassInvariant;
use vectorspace::VectorSpace;
use tol::Param;
use itertools::Itertools;
use curve::Curve;
use smat::Smat;
pub struct Bezier<P: VectorSpace> {
    spl: Bspline<P>,
}


impl<P: VectorSpace> ClassInvariant for Bezier<P> {
    fn is_valid(&self) -> Result<bool, &str> {
        try!(self.to_spline().is_valid());
        let sz = self.degree() as usize + 1;
        let t = self.knots();
        if t.len() != 2 * sz || self.control_points().len() != sz {
            return Err("Bezier knots are invalid");
        }

        if self.start_mult() == sz && self.end_mult() == sz {
            Ok(true)
        } else {
            Err("Bezier knots are invalid")
        }
    }
}


impl<P: VectorSpace> SplineWrapper for Bezier<P> {
    type TW = P;
    fn to_spline(&self) -> &Bspline<Self::TW> {
        &self.spl
    }

    fn from_spline(spl: Bspline<P>) -> Self {
        let bz = Bezier { spl: spl };
        assert!(bz.is_valid().is_ok());
        bz
    }
}

impl<P: VectorSpace> SplineMut for Bezier<P> {
    fn into_spline(self) -> Bspline<Self::T> {
        self.spl
    }
}

impl<P: VectorSpace> Curve for Bezier<P> {
    type T = P;
    fn param_range(&self) -> (f64, f64) {
        self.to_spline().param_range()
    }

    fn eval(&self, v: f64) -> Self::T {
        self.to_spline().eval(v)
    }

    fn eval_derivative(&self, v: f64, order: u32) -> Self::T {
        self.to_spline().eval_derivative(v, order)
    }
}



pub fn split_into_bezier_patches<SplineType>(spl: &SplineType) -> Vec<Bezier<SplineType::T>>
    where SplineType: SplineData + SplineMut
{
    let t = spl.knots();
    let tlen = spl.knots().len();
    let d = spl.degree() as usize;

    let mut uniq_ts = t[d..(tlen - d)]
                          .iter()
                          .map(|&x| -> Param { Param(x) })
                          .dedup()
                          .peekable();
    let mut patches: Vec<Bezier<SplineType::T>> = Vec::with_capacity(tlen - 2 * d);
    let cpts = spl.control_points();
    while let Some(t) = uniq_ts.next() {
        let Param(a) = t;
        if let Some(t) = uniq_ts.peek() {
            let &Param(b) = t;
            let nu = spl.locate_nu(a);
            let sm = Smat::new(a, b, spl.knots(), nu, spl.degree());
            let bzcpts = sm.seval(&cpts[nu - d..]);
            let mut ks = vec![a;d+1];
            ks.extend(vec![b;d+1]);
            patches.push(Bezier::new(bzcpts, ks));
        }
    }
    patches
}
#[test]
fn it_works() {
    use point::Point2;
    use tol::Tol;
    let pts = vec![Point2::new(0.0, 0.0),
                   Point2::new(1.0, 1.0),
                   Point2::new(1.5, 0.3),
                   Point2::new(1.8, 0.1),
                   Point2::new(2.0, 0.0)];
    let ks = vec![0.0, 0.0, 0.0, 0.5, 0.8, 1.0, 1.0, 1.0];
    let bs = Bspline::new(pts, ks);
    let us = vec![0.1, 0.6, 0.9];
    let mut j: usize = 0;
    for bp in split_into_bezier_patches(&bs).iter() {
        assert!(bp.is_valid().is_ok());
        assert!((bp.eval(us[j]) - bs.eval(us[j])).len().small());
        j += 1;
    }
  
}
