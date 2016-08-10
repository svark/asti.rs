use std::ops::{Add, Sub, Mul};
pub trait VectorSpace : Mul<f64, Output = Self> + Copy + Add<Self, Output=Self> + Sub<Self, Output=Self> + Default
{
    fn lerp(&self, l: f64, q: Self) -> Self {
        let p = *self;
        let m = 1.0 - l;
        p * m + q * l
    }

    fn dlerp(&self, l: f64, q: Self) -> Self {
        let p = *self;
        (q - p) * l
    }

    fn dim() -> u32;
    fn load(x: &[f64], i: usize) -> Self;
    fn splat(x: f64) -> Self;
    fn extract(&self, i: u32) -> f64;
    fn replace(&mut self, i: u32, x: f64) -> Self;

    type L : VectorSpace;
    type H : VectorSpace;
    fn ldim(&self) -> Self::L;
    fn hdim(&self, pad: f64) -> Self::H;
}
