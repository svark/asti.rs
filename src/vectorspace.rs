use std::marker::Copy;
use std::clone::Clone;
use tol::{Tol, RESABS};
use std::fmt::Debug;
use std::ops::{Add, Sub, Mul, Div,  AddAssign, MulAssign, DivAssign, Index};
use nalgebra::{Repeat, NumVector, Indexable, Dimension, Dot, Norm, Origin};
use std::mem;

pub use nalgebra::Cross;
pub use nalgebra::Axpy;

pub trait Ops : Repeat<f64> + Indexable<usize,f64> + Dimension {
    fn splat(x: f64) -> Self {
        Self::repeat(x)
    }
    fn dim() -> usize {
        Self::dimension(None)
    }

    fn extract(&self, i: usize) -> f64 {
        debug_assert!(i < Self::dim());
        unsafe { self.unsafe_at(i) }
    }
    fn replace(&mut self, i: usize, x: f64) {
        debug_assert!(i < Self::dim());
        unsafe { self.unsafe_set(i, x) }
    }
}

pub trait NVS : Norm<NormType=f64> + Dot<f64> + Mul<f64,Output=Self>  {
    fn norm(&self) -> f64 {
        Norm::norm(self)
    }
    
    fn dot(&self, other: &Self) -> f64 {
        Dot::dot(self,other)
    }

    fn norm_squared(&self) -> f64 {
        Norm::norm_squared(self)
    }
    
    fn try_normalize(&self) -> Option<Self> {
        Norm::try_normalize(self, RESABS)        
    }
}

pub trait PV : Indexable<usize,f64> + Sized + Clone {
   type V: NumVector<f64> + NVS + Indexable<usize,f64> + Clone + Copy;

   fn to_vector(&self) -> Self::V {
        let v: &Self::V = unsafe {
            mem::transmute(self)
        };
        v.clone()
   }

   fn as_vector(&self) -> &Self::V {
        unsafe {
            mem::transmute(self)
        }
   }
   
   fn from_vec(v:&Self::V) -> Self
   {
       let p: &Self = unsafe {
           mem::transmute(v)
       };
       p.clone()
   }
}


pub trait PointT :  PV + Ops  + Sized +  Copy + Clone + Debug 
  + Dimension + Origin + PartialEq + Axpy<f64> 
  + Sub<Self, Output=<Self as PV>::V> + Mul<f64, Output=Self> + Div<f64, Output=Self> + Add< <Self as PV>::V, Output=Self> + MulAssign<f64> + DivAssign<f64> 
  + AddAssign<<Self as PV>::V> + Index<usize, Output=f64>
{
    type L : PointT;
    type H : PointT;

    fn ldim(&self) -> Self::L;

    fn proj(&self) -> Option<Self> {
        let v = self.extract(Self::dim() - 1);
        if !v.small() {
            Some((*self) * (1.0 / v))
        } else {
            None
        }
    }

    fn hdim(&self, pad: f64) -> Self::H;

    fn lerp(&self, l: f64, q: Self) -> Self;

    fn zero_pt() -> Self {
        Self::splat(0.0)
    }

}



pub fn to_pt<P: PointT>(v: <P as PV>::V) -> P
{
    P::zero_pt() + v
}

pub trait AssocPoint
{
   type TW : PointT;
}

macro_rules! TW {
    () =>  {<Self as AssocPoint>::TW}
}

macro_rules! TWL {
    () =>  {<< Self as AssocPoint>::TW as PointT>::L}
}
