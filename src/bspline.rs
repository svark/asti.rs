use tol::{Tol, PARAMRES};
use vectorspace::PointT;
use util::merge;
use splinedata::{SplineData, KnotManip};
use curve::{Curve, FiniteCurve, BlossomCurve};
use rmat::{eval, locate_nu};
use nalgebra::Axpy;
use class_invariant::ClassInvariant;
use std::fmt;

#[derive(Clone,Debug)]
pub struct Bspline<Point: PointT> {
    control_points: Vec<Point>,
    knots: Vec<f64>,
    deg: u32,
}


impl<P> fmt::Display for Bspline<P>
    where P: PointT
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        try!(writeln!(f, "cpts:"));
        for &c in self.control_points() {
            try!(write!(f, "("));
            for j in 0..P::dim() {
                try!(write!(f, "{:?},", c.extract(j)));
            }
            try!(writeln!(f, ")"));
        }
        try!(write!(f, "["));
        for t in self.knots() {
            try!(write!(f, "{}", t));
        }
        write!(f, "{}", self.degree())
    }
}

pub trait SplineWrapper
{
    type TW:PointT;
    fn to_spline(&self) -> &Bspline<Self::TW>;
    fn from_spline(spl: Bspline<Self::TW>) -> Self;
}

impl<P: PointT> SplineWrapper for Bspline<P> {
    type TW = P;
    fn to_spline(&self) -> &Bspline<P> {
        self
    }
    fn from_spline(spl: Bspline<Self::TW>) -> Self {
        spl
    }
}

macro_rules! TW {
    () =>  {<Self as SplineWrapper>::TW}
}

macro_rules! TWL {
    () =>  {<< Self as SplineWrapper>::TW as PointT>::L}
}

impl<SplType> KnotManip for SplType
    where SplType: SplineWrapper
{
    fn start_mult(&self) -> usize {
        let f = self.front();
        self.knots()
            .iter()
            .take_while(|&x| (*x - f).small_param())
            .count()
    }

    fn end_mult(&self) -> usize {
        let l = self.back();
        self.knots()
            .iter()
            .rev()
            .take_while(|&x| (*x - l).small_param())
            .count()
    }

    fn front(&self) -> f64 {
        self.knots()[self.degree() as usize]
    }

    fn back(&self) -> f64 {
        self.knots()[self.knots().len() - 1 - (self.degree() as usize)]
    }

    fn locate_nu(&self, u: f64) -> usize {
        locate_nu(u, self.degree() as usize, self.knots())
    }

    fn mult(&self, u: f64) -> usize {
        let f = self.front();
        let b = self.back();
        if (u - b).small() {
            self.end_mult()
        } else if (u - f).small() {
            self.start_mult()
        } else if u > b {
            0
        } else if u < f {
            0
        } else {
            let nu = self.locate_nu(u);
            self.knots()[..nu + 1]
                .iter()
                .rev()
                .take_while(|&x| (*x - u).small_param())
                .count()
        }
    }

    fn rebase(&self, taus: Vec<f64>) -> Self {
        let d = self.degree() as usize;
        let ncpts = taus.len() - d - 1;
        let mut cpts: Vec<<Self as SplineData>::T> = Vec::with_capacity(ncpts);
        for i in 0..ncpts {
            cpts.push(self.blossom_eval(0, &taus[i..]));
        }
        Self::from_spline(Bspline::new(cpts, taus))
    }

    fn insert_knot(&self, tau: f64) -> Self {
        let nu = self.locate_nu(tau);
        let mut newts: Vec<f64> = self.knots()[0..nu + 1].to_owned();
        newts.push(tau);
        newts.extend(self.knots()[nu + 1..].iter());
        self.rebase(newts)
    }

    fn insert_knots(&self, taus: &Vec<f64>) -> Self {
        let merged = merge(self.knots(), taus);
        self.rebase(merged)
    }
}

impl<SplType> BlossomCurve for SplType
    where SplType: SplineData
{
    type T = <Self as SplineData>::T;
    fn blossom_eval(&self, der_order: u32, us: &[f64]) -> Self::T {
        let d = self.degree() as usize;
        let nu = locate_nu(us[0], d, self.knots());
        let basis = eval(self.knots(), (us, Some(nu)), self.degree(), der_order);
        let mut r: Self::T = Self::T::zero_pt();
        let cpts = &self.control_points()[nu - d..];
        for (x, y) in cpts.iter().zip(basis.iter()) {
            r.axpy(y, x);
        }
        r
    }
}


impl<P: PointT> ClassInvariant for Bspline<P> {
    fn is_valid(&self) -> Result<bool, &str> {
        let d = self.degree() as usize;
        let cptslen = self.control_points().len();
        if self.knots().len() != cptslen + d + 1 {
            return Err("degree != #knots - #cpts - 1 ");
        }
        if self.degree() < 1 {
            return Err("zero degree");
        }

        if self.control_points().len() < 1 {
            return Err("too few cpts");
        }

        let t = self.knots();
        for j in 1..t.len() {
            let i = j - 1;
            if t[i] > t[j] + PARAMRES / 2.0 {
                return Err("knots not sorted");
            }
        }
        let mut nextj: usize = 0;
        for j in 0..t.len() {
            if j < nextj {
                continue;
            }

            let m = self.mult(self.knots[j]);
            if m > self.degree() as usize + 1 {
                return Err("too many dups in knots");
            }

            nextj = j + m;
        }
        return Ok(true);
    }
}


impl<SplType: SplineWrapper> SplineData for SplType {
    type T = <Self as SplineWrapper>::TW;
    fn control_points(&self) -> &Vec<Self::T> {
        &self.to_spline().control_points
    }

    fn knots(&self) -> &Vec<f64> {
        &self.to_spline().knots
    }

    fn degree(&self) -> u32 {
        self.to_spline().deg
    }
}


impl<P> Curve for Bspline<P>
    where P: PointT
{
    type T = P;

    fn eval(&self, u: f64) -> P {
        let d = self.deg as usize;
        self.blossom_eval(0, &vec![u;d+1])
    }

    fn eval_derivative(&self, u: f64, order: u32) -> P {
        let d = self.deg as usize;
        self.blossom_eval(order, &vec![u;d+1])
    }
}

impl<P: PointT> FiniteCurve for Bspline<P> {
    fn param_range(&self) -> (f64, f64) {
        let t = &self.knots;
        let d = self.deg as usize;
        let ncpts = t.len() - d - 1;
        (t[d], t[ncpts])
    }
}
pub trait SplineMut: SplineData
{
    fn into_spline(self) -> Bspline<Self::T>;
}


impl<P: PointT> SplineMut for Bspline<P> {
    fn into_spline(self) -> Bspline<Self::T> {
        self
    }
}


impl<P: PointT> Bspline<P> {
    pub fn new(cpts: Vec<P>, ks: Vec<f64>) -> Bspline<P> {
        let d = ks.len() - cpts.len() - 1;
        Bspline {
            control_points: cpts,
            knots: ks,
            deg: d as u32,
        }
    }
}
// revisit internal api
pub fn change_knots<P: PointT>(ks: Vec<f64>, spl: Bspline<P>) -> Bspline<P> {
    debug_assert!({
        let d = spl.degree() as usize;
        let kstart = ks[d];
        let scale = (ks[ks.len() - d - 1] - kstart) / (spl.end_param() - spl.start_param());
        let mut not_matching = false;
        for (t, u) in ks.iter().zip(spl.knots().iter()) {
            if !((t - kstart) - scale * (u - spl.start_param())).small_param() {
                not_matching = true;
                break;
            }
        }
        !not_matching
    });
    let s = Bspline { knots: ks, ..spl };
    s
}

#[test]
fn it_works() {
    use point::Pt1;
    let bs = Bspline {
        control_points: vec![Pt1::new(0.0), Pt1::new(1.0), Pt1::new(0.5)],
        knots: vec![1.0, 1.0, 1.0, 2.0, 2.0, 2.0],
        deg: 2,
    };
    assert_eq!(bs.is_valid(), Ok(true));
    let vs = vec![Pt1::new(0.0), Pt1::new(0.625), Pt1::new(0.5)];
    assert_eq!(bs.eval(1.0), vs[0]);
    assert_eq!(bs.eval(1.5), vs[1]);
    assert_eq!(bs.eval(2.0), vs[2]);

    assert_eq!(bs.blossom_eval(0, &[1.0, 1.0, 1.0]), Pt1::new(0.0));
    assert_eq!(bs.blossom_eval(0, &[2.0, 2.0, 2.0]), Pt1::new(0.5));
    assert_eq!(bs.blossom_eval(0, &[1.0, 1.0, 2.0]), Pt1::new(1.0));

    assert_eq!(bs.eval_derivative(1.0, 1), Pt1::new(2.0));

    let spl = bs.insert_knot(1.1);
    assert_eq!(spl.is_valid(), Ok(true));
    assert_eq!(spl.eval(1.0), vs[0]);
    assert_eq!(spl.eval(1.5), vs[1]);
    assert_eq!(spl.eval(2.0), vs[2]);
}
