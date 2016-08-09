use tol::PARAMRES;
use tol::Tol;

use util::{merge, upper_bound};

pub struct Bspline<Point> {
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

trait KnotManip
{
    fn start_mult(&self) -> usize;
    fn end_mult(&self) -> usize;
    fn mult(&self, u: f64) -> usize;
    fn front(&self) -> f64;
    fn back(&self) -> f64;
    fn knots(&self) -> &Vec<f64>;
    fn locate_nu(&self, u: f64) -> usize;
    fn param_range(&self) -> (f64, f64);
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

impl<Point> KnotManip for Bspline<Point> {
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

    fn knots(&self) -> &Vec<f64> {
        &self.knots
    }

    fn param_range(&self) -> (f64, f64) {
        let t = &self.knots;
        let d = self.deg as usize;
        let ncpts = self.knots.len() - d - 1;
        (t[d], t[ncpts])
    }

    fn locate_nu(&self, u: f64) -> usize {
        locate_nu(u, self.deg as usize, &self.knots)
    }


    fn mult(&self, u: f64) -> usize {
        let nu = self.locate_nu(u);
        self.knots[..nu + 1].iter().rev().take_while(|&x| (*x - u).small_param()).count()
    }
}

use std::ops::{Add, SubAssign, AddAssign, Mul};


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


impl<T> Bspline<T>
    where T: Copy + Add<T> + Mul<f64, Output = T> + AddAssign<T> + SubAssign<T> + Default
{
    pub fn new(cpts: Vec<T>, ks: Vec<f64>) -> Bspline<T> {
        let d = ks.len() - cpts.len() - 1;
        Bspline {
            control_points: cpts,
            knots: ks,
            deg: d as u32,
        }
    }

    pub fn rebase(&self, taus: Vec<f64>) -> Bspline<T> {
        let d = self.deg as usize;
        let mut cpts: Vec<T> = Vec::new();
        let ncpts = taus.len() - d - 1;
        cpts.reserve(ncpts);
        for i in 0..ncpts {
            cpts.push(self.blossom_eval(0, &taus[i..]));
        }
        Bspline::new(cpts, taus)
    }

    pub fn insert_knot(&self, tau: f64) -> Bspline<T> {
        let nu = self.locate_nu(tau);
        let mut newts: Vec<f64> = self.knots[0..nu + 1].iter().cloned().collect();
        newts.push(tau);
        newts.extend(self.knots[nu + 1..].iter());
        let spl = self.rebase(newts);
        spl
    }


    pub fn insert_knots(&self, taus: &Vec<f64>) -> Bspline<T> {
        let merged = merge(&self.knots, taus);
        self.rebase(merged)
    }

    pub fn is_valid(&self) -> Result<bool, &str> {
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

    pub fn eval(&self, u: f64) -> T {
        let d = self.deg as usize;
        self.blossom_eval(0, &vec![u;d+1])
    }

    pub fn blossom_eval(&self, der_order: u32, us: &[f64]) -> T {
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
        let mut r: T = Default::default();
        let cpts = &self.control_points[nu - d..];
        for (&x, &y) in cpts.iter().zip(b.iter()) {
            r += x * y
        }
        r
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
