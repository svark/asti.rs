use tol::PARAMRES;
use tol::Tol;
use vectorspace::VectorSpace;
use util::{merge, upper_bound};

pub struct Bspline<Point:VectorSpace> {
    control_points: Vec<Point>,
    knots: Vec<f64>,
    deg: u32,
}

#[inline]
pub fn sdiv(x: f64, y: f64) -> f64 {
    if y.small() {
        0.0f64
    } else {
        x / y
    }
}

pub trait KnotManip
{
    fn start_mult(&self) -> usize;
    fn end_mult(&self) -> usize;
    fn mult(&self, u: f64) -> usize;
    fn front(&self) -> f64;
    fn back(&self) -> f64;
    fn locate_nu(&self, u: f64) -> usize;
    fn rebase(&self, taus: Vec<f64>) -> Self;
    fn insert_knot(&self, tau: f64) -> Self;
    fn insert_knots(&self, taus: &Vec<f64>) ->Self;
}

pub trait SplineData
{
    type T : VectorSpace;
    fn control_points(&self) -> &Vec<Self::T>;
    fn knots(&self) -> &Vec<f64> ;
    fn degree(&self) -> u32;
}

impl<P> SplineData for Bspline<P> where P:VectorSpace
{
    type T = P;
    fn control_points(&self) -> &Vec<P> {
        &self.control_points
    }
    fn knots(&self) -> &Vec<f64> {
        &self.knots
    }
    fn degree(&self) -> u32 {
        self.deg
    }
}

pub fn locate_nu(u: f64, d: usize, t: &[f64]) -> usize {
    let ncpts = t.len() - d - 1;
    let mut u = u;
    if u <= t[d] {
        u = t[d] + PARAMRES / 2.0;
    }
    if u >= t[ncpts] {
        u = t[ncpts] - PARAMRES / 2.0;
    }
    let idx = upper_bound(t, u) - 1;
    assert!(idx >= d && idx < ncpts);
    idx
}

impl<T> KnotManip for Bspline<T> where T:VectorSpace {
    fn start_mult(&self) -> usize {
        let f = self.front();
        self.knots.iter().take_while(|&x| (*x - f).small_param()).count()
    }

    fn end_mult(&self) -> usize {
        let l = self.back();
        self.knots.iter().rev().take_while(|&x| (*x - l).small_param()).count()
    }

    fn front(&self) -> f64 {
        self.knots[self.deg as usize]
    }

    fn back(&self) -> f64 {
        self.knots[self.knots.len() - 1 - (self.deg as usize)]
    }


    fn locate_nu(&self, u: f64) -> usize {
        locate_nu(u, self.deg as usize, &self.knots)
    }

    fn mult(&self, u: f64) -> usize {
        let nu = self.locate_nu(u);
        self.knots[..nu + 1].iter().rev().take_while(|&x| (*x - u).small_param()).count()
    }

    fn rebase(&self, taus: Vec<f64>) -> Bspline<T> {
        let d = self.deg as usize;
        let ncpts = taus.len() - d - 1;
        let mut cpts: Vec<T> = Vec::with_capacity(ncpts);
        for i in 0..ncpts {
            cpts.push(self.blossom_eval(0, &taus[i..]));
        }
        Bspline::new(cpts, taus)
    }

    fn insert_knot(&self, tau: f64) -> Bspline<T> {
        let nu = self.locate_nu(tau);
        let mut newts: Vec<f64> = self.knots[0..nu + 1].to_owned();
        newts.push(tau);
        newts.extend(self.knots[nu + 1..].iter());
        let spl = self.rebase(newts);
        spl
    }

    fn insert_knots(&self, taus: &Vec<f64>) -> Bspline<T> {
        let merged = merge(&self.knots, taus);
        self.rebase(merged)
    }
}


struct RMatTau<'a, 'b> {
    knots: &'a [f64],
    taus: &'b [f64],
    offset: usize,
}

struct RMatExplicitDer<'a, 'b>(RMatTau<'a, 'b>);


trait Entries
{
    fn get_diag_from_ndiag(v: f64) -> f64;
    fn get_ndiag(&self, sz: usize, d: usize) -> f64;
    fn get_knots(&self) -> &[f64];

}

impl<'a, 'b> Entries for RMatTau<'a, 'b> {
    #[inline]
    fn get_ndiag(&self, sz: usize, d: usize) -> f64 {
        let t = self.knots;
        let idx = sz + self.offset;
        assert!(idx >= sz);
        let denom = t[idx] - t[idx - d];
        sdiv(self.taus[d] - t[idx - d], denom)
    }

    #[inline]
    fn get_diag_from_ndiag(v: f64) -> f64 {
        1.0 - v
    }

    #[inline]
    fn get_knots(&self) -> &[f64] {
        self.knots
    }
}


impl<'a, 'b> Entries for RMatExplicitDer<'a, 'b> {
    #[inline]
    fn get_ndiag(&self, sz: usize, d: usize) -> f64 {
        let t = self.0.knots;
        let idx = sz + self.0.offset;
        assert!(idx >= sz);
        let denom = t[idx] - t[idx - d];
        1.0 / denom
    }

    #[inline]
    fn get_diag_from_ndiag(v: f64) -> f64 {
        -v
    }

    #[inline]
    fn get_knots(&self) -> &[f64] {
        self.0.knots
    }
}

struct RMatExplicitResult {
    basis: Vec<f64>,
}

impl RMatExplicitResult {
    fn mult<T>(&mut self, mat: &T, fac: f64)
        where T: Entries
    {
        let basis = &mut self.basis;
        let num_cols = basis.len();
        let mut c = mat.get_ndiag(num_cols, num_cols);
        let e = basis[num_cols - 1] * c;
        basis.resize(num_cols + 1, e);
        for l in 1..num_cols {
            let k = num_cols - l;
            let b = T::get_diag_from_ndiag(c);
            c = mat.get_ndiag(k, num_cols);
            basis[k] = fac * (b * basis[k] + c * basis[k - 1]);
        }
        basis[0] *= fac * T::get_diag_from_ndiag(c);
    }
    fn drain(self) -> Vec<f64> {
        self.basis
    }
}

pub trait ClassInvariant
{
    fn is_valid(&self) -> Result<bool,&str>;
}

impl<P> ClassInvariant for Bspline<P> where P:VectorSpace
{
    fn is_valid(&self) -> Result<bool, &str> {
        if self.knots.len() != self.control_points.len() + self.deg as usize + 1 {
            return Err("bad degree");
        }
        if self.deg < 1 {

            return Err("zero degree");
        }

        if self.control_points.len() < 1 {
            return Err("too few cpts");
        }

        let t = &self.knots;
        for j in 1..t.len() {
            let i = j - 1;
            if t[i] > t[j] {
                return Err("knots not sorted");
            }
        }
        let mut nextj: usize = 0;
        for j in 0..t.len() {
            if j < nextj {
                continue;
            }

            let m = self.mult(self.knots[j]);
            if m > self.deg as usize + 1 {
                return Err("too many dups in knots");
            }

            nextj = j + m;
        }
        return Ok(true);
    }

}
pub trait Curve
{
    type T:VectorSpace;
    fn param_range(&self)->(f64,f64);
    fn eval(&self,v:f64) -> Self::T;
    fn eval_derivative(&self, v:f64,order:u32) -> Self::T;
}

pub trait BlossomCurve
{
    type T : VectorSpace;
    fn blossom_eval(&self, der_order: u32, us: &[f64]) -> Self::T;
}

impl<P> BlossomCurve for Bspline<P> where P:VectorSpace
{
    type T = P;
    fn blossom_eval(&self, der_order: u32, us: &[f64]) -> P {
        let nu = self.locate_nu(us[0]);
        let d = self.deg as usize;
        let rmat = RMatTau {
            knots: &self.knots.as_slice(),
            taus: us,
            offset: nu,
        };
        let b = {
            let mut res = RMatExplicitResult { basis: vec![1.0f64] };
            res.basis.reserve(self.deg as usize);

            for _ in 1..(self.deg - der_order + 1) {
                res.mult(&rmat, 1.0);
            }

            let rmat_der = RMatExplicitDer(rmat);
            for j in (self.deg - der_order + 1)..(self.deg + 1) {
                res.mult(&rmat_der, j as f64);
            }
            res.drain()
        };
        let mut r: P = Default::default();
        let cpts = &self.control_points[nu - d..];
        for (&x, &y) in cpts.iter().zip(b.iter()) {
            r = r + x * y
        }
        r
    }
}


impl<P> Curve for Bspline<P> where P:VectorSpace
{
    type T = P;
    fn param_range(&self) -> (f64, f64) {
        let t = &self.knots;
        let d = self.deg as usize;
        let ncpts = t.len() - d - 1;
        (t[d], t[ncpts])
    }

    fn eval(&self, u: f64) -> P {
        let d = self.deg as usize;
        self.blossom_eval(0, &vec![u;d+1])
    }

    fn eval_derivative(&self, u: f64, order: u32) -> P
    {
        let d = self.deg as usize;
        self.blossom_eval(order, &vec![u;d+1])
    }
}

impl<T> Bspline<T>
    where T: VectorSpace
{
    pub fn new(cpts: Vec<T>, ks: Vec<f64>) -> Bspline<T> {
        let d = ks.len() - cpts.len() - 1;
        Bspline {
            control_points: cpts,
            knots: ks,
            deg: d as u32,
        }
    }
}

#[test]
fn it_works() {
    let bs = Bspline {
        control_points: vec![0.0, 1.0, 0.5],
        knots: vec![1.0, 1.0, 1.0, 2.0, 2.0, 2.0],
        deg: 2,
    };
    assert_eq!(bs.is_valid(), Ok(true));
    assert_eq!(bs.eval(1.0), 0.0);
    assert_eq!(bs.eval(1.5), 0.625);
    assert_eq!(bs.eval(2.0), 0.5);

    assert_eq!(bs.blossom_eval(0, &[1.0, 1.0, 1.0]), 0.0);
    assert_eq!(bs.blossom_eval(0, &[2.0, 2.0, 2.0]), 0.5);
    assert_eq!(bs.blossom_eval(0, &[1.0, 1.0, 2.0]), 1.0);

    let spl = bs.insert_knot(1.1);
    assert_eq!(spl.is_valid(), Ok(true));
    assert_eq!(spl.eval(1.0), 0.0);
    assert_eq!(spl.eval(1.5), 0.625);
    assert_eq!(spl.eval(2.0), 0.5);
}
