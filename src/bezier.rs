use splinedata::{SplineData, KnotManip};
use bspline::{Bspline, SplineWrapper, SplineMut};
use class_invariant::ClassInvariant;
use vectorspace::PointT;
use tol::Param;
use itertools::Itertools;
use curve::{Curve, FiniteCurve};
use smat::Smat;

#[derive(Debug)]
pub struct Bezier<P: PointT> {
    spl: Bspline<P>,
}

impl<P: PointT> Bezier<P> {
    pub fn new(cpts: Vec<P>, s: f64, e: f64) -> Bezier<P> {
        let mut knots = vec![s;cpts.len()];
        knots.extend(vec![e;cpts.len()]);
        Bezier { spl: Bspline::new(cpts, knots) }
    }
}

impl<P: PointT> ClassInvariant for Bezier<P> {
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


impl<P: PointT> SplineWrapper for Bezier<P> {
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

impl<P: PointT> SplineMut for Bezier<P> {
    fn into_spline(self) -> Bspline<Self::T> {
        self.spl
    }
}

impl<P: PointT> Curve for Bezier<P> {
    type T = P;

    fn eval(&self, v: f64) -> Self::T {
        self.to_spline().eval(v)
    }

    fn eval_derivative(&self, v: f64, order: u32) -> Self::T {
        self.to_spline().eval_derivative(v, order)
    }
}

impl<P: PointT> FiniteCurve for Bezier<P> {
    fn param_range(&self) -> (f64, f64) {
        self.to_spline().param_range()
    }
}


pub fn split_into_bezier_patches<SplineType>(spl: &SplineType) -> Vec<Bezier<SplineType::T>>
    where SplineType: SplineWrapper + SplineMut
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
            patches.push(Bezier::from_spline(Bspline::new(bzcpts, ks)));
        }
    }
    patches
}
#[test]
fn it_works() {
    use point::Pt2;
    use tol::Tol;
    use nalgebra::Norm;
    let pts = vec![Pt2::new(0.0, 0.0),
                   Pt2::new(1.0, 1.0),
                   Pt2::new(1.5, 0.3),
                   Pt2::new(1.8, 0.1),
                   Pt2::new(2.0, 0.0)];
    let ks = vec![0.0, 0.0, 0.0, 0.5, 0.8, 1.0, 1.0, 1.0];
    let bs = Bspline::new(pts, ks);
    let us = vec![0.1, 0.6, 0.9];
    let mut j: usize = 0;
    for bp in split_into_bezier_patches(&bs).iter() {
        assert!(bp.is_valid().is_ok());
        assert!((bp.eval(us[j]) - bs.eval(us[j])).norm().small());
        j += 1;
    }

}
