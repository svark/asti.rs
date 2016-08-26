use std::ops::{Add, Sub, Mul};
use std::marker::Copy;
use std::clone::Clone;
use tol::Tol;
pub trait VectorSpace : Mul<f64, Output = Self> + Copy + Add<Self, Output=Self> + Sub<Self, Output=Self> + Default + Clone
{
    fn lerp(&self, l: f64, q: Self) -> Self {
        let p = *self;
        let m = 1.0 - l;
        p * m + q * l
    }

    fn dim() -> u32;
    fn load(x: &[f64], i: usize) -> Self;
    fn splat(x: f64) -> Self;
    fn extract(&self, i: u32) -> f64;
    fn replace(&mut self, i: u32, x: f64) -> Self;

    type L : VectorSpace;
    type H : VectorSpace;
    fn ldim(&self) -> Self::L;
    fn proj(&self) -> Option<Self> {
        let v = self.extract(Self::dim());
        if !v.small() {
            Some((*self) *(1.0/v))
        }else {
            None
        }
    }
    fn hdim(&self, pad: f64) -> Self::H;
    fn len(&self) -> f64 {  self.dot(self).sqrt()  }
    fn lensq(&self) -> f64 {  self.dot(self)  }
    fn dot(&self, v: &Self) -> f64;

    fn normalize(&self) -> Option<Self> {
        if self.len().small()
        {
            None
        }else {
            Some(*self * (1.0/self.len()))
        }
    }
}

