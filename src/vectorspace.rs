use std::marker::Copy;
use std::clone::Clone;
use tol::Tol;
use std::fmt::Debug;
use std::ops::Add;
use nalgebra::{Repeat, NumPoint, Indexable, Dimension, PointAsVector};

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


pub trait PointT :   Ops +  NumPoint<f64> + Sized +  Copy + Clone + Debug
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

pub fn to_pt<P: PointT>(v: <P as PointAsVector>::Vector) -> P
    where P: Add<<P as PointAsVector>::Vector>
{
    P::zero_pt() + v
}
