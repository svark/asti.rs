use std::marker::Copy;
use std::clone::Clone;
use tol::Tol;
use std::fmt::Debug;
use nalgebra::{Repeat, NumPoint, Indexable, Dimension, Norm, Dot, FloatVector, Origin};

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



pub trait PointT :   Ops +  NumPoint<f64> + Sized + Origin + Copy + Clone + Debug
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

pub trait VectorT :   Ops + Norm<NormType=f64>  + Dot<f64> + FloatVector<f64>  {
    type L : VectorT;
    type H : VectorT;
    type P : PointT;

    fn ldim(&self) -> Self::L;
    fn hdim(&self, pad: f64) -> Self::H;

    fn to_pt(&self) -> Self::P;
}
