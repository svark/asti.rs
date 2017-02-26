use splinedata::{SplineData, KnotManip};
use bspline::{Bspline, SplineWrapper};
use curve::{Curve, Domain};
use vectorspace::{PointT, AssocPoint, NVS};
use tol::Tol;
use class_invariant::ClassInvariant;
pub struct PeriodicBspline<Point: PointT> {
    spl: Bspline<Point>,
}


impl<P: PointT> AsRef<Bspline<P>> for PeriodicBspline<P> {
    fn as_ref(&self) -> &Bspline<P> {
        &self.spl
    }
}

impl<P: PointT> AsMut<Bspline<P>> for PeriodicBspline<P> {
    fn as_mut(&mut self) -> &mut Bspline<P> {
        &mut self.spl
    }
}

impl<P: PointT> From<Bspline<P>> for PeriodicBspline<P> {
    fn from(spl: Bspline<P>) -> PeriodicBspline<P> {
        PeriodicBspline { spl: spl }
    }
}

impl<P: PointT> AssocPoint for PeriodicBspline<P> {
    type TW = P;
}

impl<P: PointT> SplineWrapper for PeriodicBspline<P> {}

impl<Point: PointT> PeriodicBspline<Point> {
    pub fn new_from_unwrapped(pts: Vec<Point>, ks: Vec<f64>, degree: usize) -> Self {
        let mut cpts: Vec<Point> = Vec::with_capacity(pts.len() + degree);
        let sz = ks.len() - 1;
        assert!(sz == pts.len());
        cpts.extend(pts[sz - degree..sz].iter());
        cpts.extend(pts[..sz].into_iter());

        let mut t: Vec<f64> = Vec::with_capacity(sz + 1 + 2 * degree);

        let mut start_diffs = Vec::with_capacity(degree);
        let mut end_diffs = Vec::with_capacity(degree);

        for wins in ks.windows(2).take(degree) {
            start_diffs.push(wins[1] - wins[0]);
        }
        for wins in ks.windows(2).rev().take(degree) {
            end_diffs.push(wins[1] - wins[0]);
        }

        let mut f = *ks.first().unwrap();
        let mut fm: Vec<f64> = Vec::with_capacity(degree);
        for e in end_diffs.iter() {
            f -= *e;
            fm.push(f);
        }
        t.extend(fm.iter().rev());
        let mut l = *ks.last().unwrap();
        let mut bp: Vec<f64> = Vec::with_capacity(degree);
        for s in start_diffs.iter() {
            l += *s;
            bp.push(l);
        }
        t.extend(ks.iter());
        t.extend(bp.iter());
        PeriodicBspline::from(Bspline::new(cpts, t))
    }
}

pub fn periodic_param(rng: (f64, f64), u: f64) -> f64 {
    let (s, e) = rng;
    let mut r = u - s;
    let d = e - s;
    if r > 0.0 && r < d {
        u
    } else {
        r %= d;
        if r < 0.0 {
            r += d;
        }
        s + r
    }
}

impl<Point: PointT> Curve for PeriodicBspline<Point> {
    type T = Point;
    fn eval(&self, u: f64) -> Point {
        let v = periodic_param(self.param_range(), u);
        self.spl.eval(v)
    }
    fn eval_derivative(&self, u: f64, der_order: u32) -> Point {
        let v = periodic_param(self.param_range(), u);
        self.spl.eval_derivative(v, der_order)
    }
}

impl<Point: PointT> ClassInvariant for PeriodicBspline<Point> {
    fn is_valid(&self) -> Result<bool, &str> {
        try!(self.spl.is_valid());
        let pr = self.param_range();

        let mut its_periodic = true;
        let mult = self.mult(self.knots()[self.degree() as usize]);

        for i in 0..(self.degree() + 1 - (mult as u32)) {
            let v1 = self.eval_derivative(pr.0, i);
            let v2 = self.eval_derivative(pr.1, i);
            if !(v1 - v2).norm().small() {
                its_periodic = false;
                break;
            }
        }
        if its_periodic {
            Ok(true)
        } else {
            Err("Periodicity violated")
        }
    }
}

impl<P: PointT> Into<Bspline<P>> for PeriodicBspline<P> {
    fn into(self) -> Bspline<P> {
        self.spl
    }
}

#[test]
fn it_works() {
    use point::Pt2;
    let pts = vec![Pt2::new(0.0, 0.0), Pt2::new(0.4, 0.3), Pt2::new(0.2, 0.8), Pt2::new(-0.2, 0.4)];
    let ks = vec![0.0, 0.3, 0.6, 0.8, 1.0];

    let bs = PeriodicBspline::new_from_unwrapped(pts, ks, 2);
    for t in bs.knots().iter() {
        println!(">{}", t);
    }
    let (a, b) = bs.param_range();
    assert!(bs.is_valid().is_ok());
    assert_eq!(a, 0.0);
    assert_eq!(b, 1.0);
    assert!((bs.eval(0.0) - Pt2::new(0.04, 0.64)).norm().small());
    assert!((bs.eval(0.2) - Pt2::new(-0.12888889, 0.33777778)).norm().small());
    assert!((bs.eval(0.9) - Pt2::new(0.185, 0.6975)).norm().small());
}
