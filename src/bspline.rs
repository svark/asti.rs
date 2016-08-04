//extern crate tol;
use tol::PARAMRES;
use tol::Tol;


pub struct Bspline<Point> {
    control_points: Vec<Point>,
    knots: Vec<f64>,
    deg : u32
}

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

// struct RMat<'a>
// {
//     knots: &'a [f64],
//     d: u32,
// }

use std::ops::{Add, SubAssign, AddAssign, Mul};

// use point::{dlerp, lerp};
// impl<'a> RMat<'a>
// {
//     fn reval<T: Mul<f64> + Copy + Add<T> + AddAssign<T> + SubAssign<T> >
//         (&self, nu:usize, u:f64, der_order:u32, cache:&mut [T]) -> bool
//     {
//         assert!(self.knots[0] == u );
//         let t =  self.knots;
//         let mut fac : usize = 1;
//         let size = self.d;
//         if der_order > size {
//             return true;
//         }

//         for d in (size - der_order)..size {
//             let sz  = (size - d + 1) as usize;
//             fac *= sz;
//             for j in 1..(sz+1) {
//                 let b = t[j] - t[j - sz];
//                 if b.small_param() {
//                     return false;
//                 }
//                 let lambda =  1 / b;
//                 cache[j - 1] = dlerp(lambda, cache[j - 1], cache[j]);
//             }
//         }

//         for d in 0..(size - der_order) {
//             let sz = (size - d + 1) as usize;
//             for j in 1..(sz + 1) {
//                 let b = t[j] - t[j - sz];
//                 let lambda =  sdiv(u  - t[j - sz], b);
//                 cache[j - 1] = lerp(lambda, cache[j - 1], cache[j] );
//             }
//         }
//         true
//     }
// }

struct RMatExplicit<'a>
{
    knots: &'a [f64],
    d: u32,
    tau: f64,
    offset: usize
}

struct RMatExplicitDer<'a>(RMatExplicit<'a>);

// static[ 	]+\([^ ]+\)[ 	]+\([^ 	]+\)(\([^\)]*\))
// \([^ ]+\)[ 	]+\([^ 	]+\)(\([^\)]*\))

impl <'a> RMatExplicit<'a>
{
    fn size(&'a self) -> u32
    {
        return self.d + 1;
    }
}

trait Entries
{
    fn get_diag_from_ndiag(v: f64) -> f64;
    fn get_ndiag(&self, i: usize) -> f64;
    fn size(&self) -> usize;
    fn get_knots(&self) -> &[f64];

}

impl<'a> Entries for RMatExplicit<'a>
{
    #[inline]
    fn get_ndiag(&self, i: usize) -> f64
    {
        let k = self.d as usize;
        let t = self.knots;
        let idx = i + 1 + self.offset;
        assert!(idx >= k);
        let denom = t[idx] - t[idx - k];
        sdiv(self.tau - t[idx - k], denom)
    }

    #[inline]
    fn get_diag_from_ndiag(v: f64) -> f64
    {
        1.0 - v
    }

    #[inline]
    fn size(&self) -> usize { self.d as usize + 1  }

    #[inline]
    fn get_knots(&self) -> &[f64] { self.knots }
}


impl<'a> Entries for RMatExplicitDer<'a>
{
    #[inline]
    fn get_ndiag(&self, i: usize) -> f64
    {
        let k = self.0.d as usize;
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
    fn size(&self) -> usize{ self.0.d as usize + 1  }

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
        assert_eq!(mat.size(), num_cols + 1);
        let mut c = mat.get_ndiag(num_cols-1);
        let e = basis[num_cols-1] * c;
        basis.resize(num_cols + 1 , e);
        for l in 1..num_cols {
            let k = num_cols - l;
            let b = T::get_diag_from_ndiag(c);
            c = mat.get_ndiag(k-1);
            basis[k]  = b * basis[k]  + c * basis[k-1];
        }
        basis[0] *= T::get_diag_from_ndiag(c);
    }

    fn get(&self, k: usize) -> f64
    {
        self.basis[k]
    }

    fn size(&self) -> usize
    {
        self.basis.len()
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

    fn get_basis(&self, nu: usize, u: f64, der_order: u32) -> Vec<f64>
    {
        let mut res = RMatExplicitResult {
            basis: vec![1.0f64]
        };
        res.basis.reserve( self.deg as usize );
        let mut rmat_noder = RMatExplicit{
            knots : &self.knots.as_slice(),
            d: 0, tau: u, offset : nu
        };
        for j in 1..(self.deg - der_order+1)
        {
            rmat_noder.d = j;
            res.mult(&rmat_noder);
        }
        let mut rmat_der = RMatExplicitDer(rmat_noder);
        for j in (self.deg - der_order + 1)..(self.deg+1)
        {
            rmat_der.0.d = j;
            res.mult( &rmat_der);
        }
        res.basis
   }

    pub fn eval(&self, u : f64) -> T {
        let nu  = self.locate_nu(u)  ;
        let b = self.get_basis(nu, u, 0u32);
        let mut r:T =  Default::default();
        let d = self.deg as usize;
        for j in 0..d+1
        {
            let p = self.control_points[(nu - d + j) ] * b[j];
            r += p;
        }
        r
    }

    // fn modify(&self, modf : FnMut<Vec<Point>, Vec<f64> > ) -> Self
    // {
    //     (cpts, ks) = modf(self.control_points, self.knots);
    //     Bspline::new(cpts, ks);
    // }

}

#[test]
fn it_works()
{
    let bs = Bspline {control_points: vec![0.0,1.0,0.5], knots : vec![1.0,1.0,1.0,2.0,2.0,2.0], deg : 2 };
    assert_eq!(bs.eval(1.0), 0.0);
    assert_eq!(bs.eval(2.0), 0.5);
}
