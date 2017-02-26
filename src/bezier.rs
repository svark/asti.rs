use splinedata::{SplineData, KnotManip};
use bspline::{Bspline, SplineWrapper};
use class_invariant::ClassInvariant;
use vectorspace::{PointT, AssocPoint};
use tol::Param;
use itertools::Itertools;
use curve::Curve;
use smat::Smat;
use std::ops::Deref;

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
        try!(self.spl.is_valid());
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

impl<P: PointT> From<Bspline<P>> for Bezier<P> {
    fn from(spl: Bspline<P>) -> Bezier<P> {
        let bz = Bezier { spl: spl };
        assert!(bz.is_valid().is_ok());
        bz
    }
}

impl<P: PointT> AsRef<Bspline<P>> for Bezier<P> {
    fn as_ref(&self) -> &Bspline<P> {
        &self.spl
    }
}

impl<P: PointT> AsMut<Bspline<P>> for Bezier<P> {
    fn as_mut(&mut self) -> &mut Bspline<P> {
        &mut self.spl
    }
}

impl<P: PointT> Into<Bspline<P>> for Bezier<P> {
    fn into(self) -> Bspline<P> {
        self.spl
    }
}

impl<P: PointT> Deref for Bezier<P> {
    type Target = Bspline<P>;
    fn deref(&self) -> &Bspline<P> {
        self.as_ref()
    }
}

impl<P: PointT> AssocPoint for Bezier<P> {
    type TW = P;
}
impl<P: PointT> SplineWrapper for Bezier<P> {}

impl<P: PointT> Curve for Bezier<P> {
    type T = P;

    fn eval(&self, v: f64) -> Self::T {
        self.as_ref().eval(v)
    }

    fn eval_derivative(&self, v: f64, order: u32) -> Self::T {
        self.as_ref().eval_derivative(v, order)
    }
}

pub fn split_into_bezier_patches<P>(spl: &Bspline<P>) -> Vec<Bezier<P>>
    where P: PointT
{
    let t = spl.knots();
    let tlen = spl.knots().len();
    let d = spl.degree() as usize;

    let mut uniq_ts = uniq_ts![t[d..(tlen - d)].iter()].peekable();
    let mut patches: Vec<Bezier<P>> = Vec::with_capacity(tlen - 2 * d);
    let cpts = spl.control_points();
    while let Some(a) = uniq_ts.next() {
        if let Some(&b) = uniq_ts.peek() {
            let nu = spl.locate_nu(a);
            let sm = Smat::new(a, b, spl.knots(), nu, spl.degree());
            let bzcpts = sm.seval(&cpts[nu - d..]);
            let mut ks = vec![a;d+1];
            ks.extend(vec![b;d+1]);
            patches.push(Bezier::from(Bspline::new(bzcpts, ks)));
        }
    }
    patches
}

#[test]
fn it_works() {
    use point::Pt2;
    use tol::Tol;
    use vectorspace::NVS;
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
