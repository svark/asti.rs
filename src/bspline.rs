use tol::PARAMRES;
use tol::Tol;

pub struct Bspline<Point> {
    control_points: Vec<Point>,
    knots: Vec<f64>,
    deg : u32
}
// its a pity rust does not have upper_bound
fn upper_bound(v: &Vec<f64>, x:f64) -> usize {
    let mut cnt:usize = v.len();
    let mut ub = 0;
    while cnt > 0 {
        let step = cnt/2;
        let mid = ub + step;
        if x >= v[mid]   {
            ub = mid + 1;
            cnt -= step + 1;
        } else {
            cnt = step;
        }
    };
    ub
}

#[inline]
fn sdiv(x: f64, y: f64) -> f64 {
    if y.small() {
        0.0f64
    }
    else {
        x / y
    }
}

trait KnotManip
{
    fn start_mult(&self) -> usize;
    fn end_mult(&self) -> usize;
    fn front(&self) -> f64;
    fn back(&self) -> f64;
    fn knots(&self) -> &Vec<f64>;
    fn locate_nu(&self,u: f64) -> usize;
}

impl <Point> KnotManip for Bspline<Point> {

    fn start_mult(&self) -> usize {
        let f = self.front();
        self.knots.iter().take_while(|&x| *x == f).count()
    }

    fn end_mult(&self) -> usize {
        let l = self.back();
        self.knots.iter().rev().take_while(|&x| *x == l).count()
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

    fn locate_nu(&self, u: f64) -> usize
    {
        let d = self.deg as usize;
        let ncpts = self.knots.len() - d - 1;
        let t = &self.knots;
        let mut u = u;
        if u <= t[d] {
            u = t[d] + PARAMRES/2.0;
        }
        if u >= t[ncpts] {
            u = t[ncpts] - PARAMRES/2.0;
        }

        let idx = upper_bound(&t,u) - 1;
        assert!(idx >= d && idx < ncpts);
        idx
    }
}

use std::ops::{Add, SubAssign, AddAssign, Mul};


struct RMatTau<'a,'b>
{
    knots: &'a [f64],
    tau: &'b [f64],
    offset: usize
}

struct RMatExplicitDer<'a,'b>(RMatTau<'a,'b>);


trait Entries
{
    fn get_diag_from_ndiag(v: f64) -> f64;
    fn get_ndiag(&self, k: usize, i: usize) -> f64;
    fn get_knots(&self) -> &[f64];

}

impl<'a,'b> Entries for RMatTau<'a,'b>
{
    #[inline]
    fn get_ndiag(&self, k: usize, i: usize) -> f64
    {
        let t = self.knots;
        let idx = i + 1 + self.offset;
        assert!(idx >= k);
        let denom = t[idx] - t[idx - k];
        sdiv(self.tau[i] - t[idx - k], denom)
    }

    #[inline]
    fn get_diag_from_ndiag(v: f64) -> f64
    {
        1.0 - v
    }

    #[inline]
    fn get_knots(&self) -> &[f64] { self.knots }
}


impl<'a,'b> Entries for RMatExplicitDer<'a,'b>
{
    #[inline]
    fn get_ndiag(&self, k: usize, i: usize) -> f64
    {
        let t = self.0.knots;
        let idx = i + 1 + self.0.offset;
        assert!(idx >= k);
        let denom = t[idx] - t[idx - k];
        1.0/denom
    }

    #[inline]
    fn get_diag_from_ndiag(v : f64) -> f64
    {
        -v
    }

    #[inline]
    fn get_knots(&self) -> &[f64] { self.0.knots }
}

struct RMatExplicitResult
{
    basis : Vec<f64>
}

impl RMatExplicitResult
{

    fn mult<T>(&mut self, mat: &T) where T: Entries
    {
        let basis = &mut self.basis;
        let num_cols = basis.len();
        let mut c = mat.get_ndiag(num_cols, num_cols-1);
        let e = basis[num_cols-1] * c;
        basis.resize(num_cols + 1 , e);
        for l in 1..num_cols {
            let k = num_cols - l;
            let b = T::get_diag_from_ndiag(c);
            c = mat.get_ndiag(k, k-1);
            basis[k]  = b * basis[k]  + c * basis[k-1];
        }
        basis[0] *= T::get_diag_from_ndiag(c);
    }
    fn drain(self) -> Vec<f64>
    {
        self.basis
    }
}


impl <T> Bspline<T> where T : Copy + Add<T> + Mul<f64, Output=T> + AddAssign<T> + SubAssign<T> + Default {
    pub fn new(cpts: Vec<T>, ks : Vec<f64>) -> Bspline<T>
    {
        let d  = ks.len() - cpts.len() -1 ;
        Bspline{ control_points : cpts, knots: ks, deg: d as u32}
    }


    fn get_blossom_basis<'a,'b>(&self, der_order : u32, rmat : RMatTau<'a,'b>) -> Vec<f64>
    {
        let mut res = RMatExplicitResult {
            basis: vec![1.0f64]
        };
        res.basis.reserve( self.deg as usize );

        for _ in 1..(self.deg - der_order+1)
        {
            res.mult(&rmat);
        }

        let rmat_der = RMatExplicitDer(rmat);
        for _ in (self.deg - der_order + 1)..(self.deg+1)
        {
            res.mult(&rmat_der);
        }
        res.drain()
    }

    pub fn eval(&self, u : f64) -> T {
        let d = self.deg as usize;
        self.blossom_eval(&vec![u;d])
    }

    pub fn blossom_eval(&self, u : &[f64]) -> T {
        let nu = self.locate_nu(u[0]);
        let d = self.deg as usize;
        let rmat = RMatTau{
            knots : &self.knots.as_slice(),
            tau: u, offset:nu
        };
        let b = self.get_blossom_basis(0, rmat);
        let mut r:T =  Default::default();
        let cpts = &self.control_points[nu-d..];
        for (&x, &y) in cpts.iter().zip(b.iter())
        {
            r += x * y
        }
        r
    }
}

#[test]
fn it_works()
{
    let bs = Bspline {control_points: vec![0.0,1.0,0.5], knots : vec![1.0,1.0,1.0,2.0,2.0,2.0], deg : 2 };
    assert_eq!(bs.eval(1.0), 0.0);
    assert_eq!(bs.eval(1.5), 0.625);
    assert_eq!(bs.eval(2.0), 0.5);

    assert_eq!(bs.blossom_eval(&[1.0,1.0]) , 0.0);
    assert_eq!(bs.blossom_eval(&[2.0,2.0]),  0.5);
//    assert_eq!(bs.blossom_eval(&[1.0,2.0]), 1.0);
}
