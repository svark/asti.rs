//extern crate tol;
use tol::PARAMRES;
use tol::Tol;
use point::{dlerp, lerp};
use std::num::Zero;

pub struct Bspline<Point> {
    control_points: Vec<Point>,
    knots: Vec<f64>,
    deg : u32
}

fn upper_bound(v: &Vec<f64>, x:f64) -> usize {
    let mut lb:usize = 0;
    let mut ub:usize = v.len();
    while lb + 1 < ub {
        let mid = (lb + ub) >> 1;
        if x < v[mid] {
            ub = mid;
        } else {
            lb = mid;
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
            u = t[ncpts] + PARAMRES/2.0;
        }

        upper_bound(&t,u) - 1
    }
}

// struct RMat<'a>
// {
//     knots: &'a [f64],
//     d: u32,
// }

use std::ops::{Add, Sub, SubAssign, AddAssign, Mul};

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
    tau: f64
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
    fn get_diag(&self, i: usize) -> f64;
    fn get_ndiag(&self, i: usize) -> f64;
    fn size(&self) -> usize;
    fn lp(lambda: f64, y: f64, z: f64) -> f64;
}

impl<'a> Entries for RMatExplicit<'a>
{
    fn get_diag(&self, i: usize) -> f64
    {
        let k = self.d as usize;
        let denom = self.knots[i + 1] - self.knots[i + 1 - k];
        sdiv(self.knots[i + 1] - self.tau, denom)
    }

    fn get_ndiag(&self, i: usize) -> f64
    {
        1.0 - self.get_diag(i)
    }

    fn size(&self) -> usize { self.d as usize + 1  }

    fn lp(lambda: f64, y: f64, z: f64) -> f64  { lerp(lambda, y, z) }
}


impl<'a> Entries for RMatExplicitDer<'a>
{
    fn get_diag(&self, i: usize) -> f64
    {
        let k = self.0.d;
        let t = self.0.knots;
        let denom = t[i + 1] - t[i + 1 - k as usize];
        -1.0/denom
    }

    fn get_ndiag(&self, i: usize) -> f64
    {
        -self.get_diag(i)
    }

    fn size(&self) -> usize{ self.0.d as usize + 1  }

    fn lp(lambda: f64, y: f64, z: f64) -> f64 { dlerp( lambda, y, z) }
}

struct RMatExplicitResult<'a>
{
    knots : &'a [f64],
    basis : Vec<f64>
}

impl<'a> RMatExplicitResult<'a>
{
    fn apply<T>(t: &Vec<f64>, mat : &T, i: usize) -> f64
        where T: Entries
    {
        let sz = t.len();
        if  i == sz
        {
            t[i-1] * mat.get_ndiag(i-1)
        }else if  i == 0
        {
            t[0] * mat.get_diag(0)
        }
        else {
            T::lp(mat.get_ndiag(i-1), t[i], t[i-1])
        }
    }

    fn mult<'b, T>(&'b mut self, kim: &T) where T: Entries
    {
        let num_cols = kim.size();
        assert_eq!(self.basis.len() + 1 , num_cols); // mult order
        let t = &mut self.basis;
        let e = t[num_cols - 2] * kim.get_ndiag(num_cols - 2);
        t.resize(num_cols, 0.0);
        for l in 0..num_cols {
            let k = num_cols - 1 - l;
            t[k] = Self::apply(t, kim, k);
        }
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
        let mut res = RMatExplicitResult { knots: &self.knots, basis : vec![1.0f64] };
        let mut rmat_noder = RMatExplicit{knots : &self.knots[nu..], d: 0, tau: u};
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
        for j in 0..d
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
