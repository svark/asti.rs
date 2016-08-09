#[cfg(target_feature = "sse2")]
use simd::x86::sse2::f64x2;

#[cfg(target_feature = "avx")]
use simd::x86::avx::f64x4;

#[cfg(target_feature = "neon")]
use simd::aarch64::neon::f64x2;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#[allow(non_camel_case_types)]
#[cfg(all(not(target_feature = "sse2"), not(target_feature = "neon")))]
#[derive(Debug, Copy, Clone)]
pub struct f64x2(f64, f64);

#[allow(non_camel_case_types)]
#[cfg(not(target_feature = "avx"))]
#[derive(Debug, Copy, Clone)]
pub struct f64x4(f64, f64, f64, f64);

#[allow(non_camel_case_types)]
#[derive(Debug, Copy, Clone)]
pub struct f64x3(f64, f64, f64);
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#[cfg(not(target_feature = "avx"))]
impl f64x4 {
    pub fn new(x: f64, y: f64, z: f64, w: f64) -> f64x4 {
        f64x4(x, y, z, w)
    }

    pub fn load(x: &[f64], i: usize) -> f64x4 {
        f64x4(x[i], x[i + 1], x[i + 2], x[i + 3])
    }

    pub fn splat(x: f64) -> f64x4 {
        f64x4(x, x, x, x)
    }

    pub fn extract(&self, i: u32) -> f64 {
        match i {
            0 => self.0,
            1 => self.1,
            2 => self.2,
            3 => self.3,
            _ => panic!("bad index for extract in f64x4"),
        }
    }

    pub fn replace(&mut self, i: u32, x: f64) -> f64x4 {
        match i {
            0 => f64x4(x, self.1, self.2, self.3),
            1 => f64x4(self.0, x, self.2, self.3),
            2 => f64x4(self.0, self.1, x, self.3),
            3 => f64x4(self.0, self.1, self.3, x),
            _ => panic!("bad index for replace in f64x4"),
        }
    }
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
impl f64x3 {
    pub fn new(x: f64, y: f64, z: f64) -> f64x3 {
        f64x3(x, y, z)
    }

    pub fn load(x: &[f64], i: usize) -> f64x3 {
        f64x3(x[i], x[i + 1], x[i + 2])
    }

    pub fn splat(x: f64) -> f64x3 {
        f64x3(x, x, x)
    }

    pub fn extract(&self, i: u32) -> f64 {
        match i {
            0 => self.0,
            1 => self.1,
            2 => self.2,
            _ => panic!("bad index for extract in f64x3"),
        }
    }

    pub fn replace(&mut self, i: u32, x: f64) -> f64x3 {
        match i {
            0 => f64x3(x, self.1, self.2),
            1 => f64x3(self.0, x, self.2),
            2 => f64x3(self.0, self.1, x),
            _ => panic!("bad index for replace in f64x3"),
        }
    }
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use std::ops::{Add, Sub, SubAssign, AddAssign, Mul};

#[cfg(not(target_feature = "avx"))]
impl Add for f64x4 {
    type Output = f64x4;

    fn add(self, p2: f64x4) -> f64x4 {
        f64x4(self.0 + p2.0, self.1 + p2.1, self.2 + p2.2, self.3 + p2.3)
    }
}

#[cfg(not(target_feature = "avx"))]
impl Sub for f64x4 {
    type Output = f64x4;

    fn sub(self, p2: f64x4) -> f64x4 {
        f64x4(self.0 - p2.0, self.1 - p2.1, self.2 - p2.2, self.3 - p2.3)
    }
}

#[cfg(not(target_feature = "avx"))]
impl AddAssign for f64x4 {
fn add_assign(&mut self, p2: f64x4) {
    *self = *self + p2;
}
}

#[cfg(not(target_feature = "avx"))]
impl SubAssign for f64x4 {
fn sub_assign(&mut self, p2: f64x4) {
    *self = *self - p2;
}
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
impl Add for f64x3 {
    type Output = f64x3;

    fn add(self, p2: f64x3) -> f64x3 {
        f64x3(self.0 + p2.0, self.1 + p2.1, self.2 + p2.2)
    }
}

impl AddAssign for f64x3 {
    fn add_assign(&mut self, p2: f64x3) {
        self.0 += p2.0;
        self.1 += p2.1;
        self.2 += p2.2;
    }
}

impl Sub for f64x3 {
    type Output = f64x3;

    fn sub(self, p2: f64x3) -> f64x3 {
        f64x3(self.0 - p2.0, self.1 - p2.1, self.2 - p2.2)
    }
}

impl SubAssign for f64x3 {
    fn sub_assign(&mut self, p2: f64x3) {
        self.0 -= p2.0;
        self.1 -= p2.1;
        self.2 -= p2.2;
    }
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#[derive(Debug, Copy, Clone)]
pub struct Point4(f64x4);

#[derive(Debug, Copy, Clone)]
pub struct Point2(f64x2);

#[derive(Debug, Copy, Clone)]
pub struct Point3(f64x3);

impl From<f64x4> for Point4 {
    fn from(v4: f64x4) -> Point4 {
        Point4(v4)
    }
}

impl From<f64x2> for Point2 {
    fn from(v2: f64x2) -> Point2 {
        Point2(v2)
    }
}

impl From<f64x3> for Point3 {
    fn from(v3: f64x3) -> Point3 {
        Point3(v3)
    }
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
impl Mul<f64> for Point4 {
    type Output = Point4;

    fn mul(self, a: f64) -> Point4 {
        let uf4 = self.0;
        Point4(f64x4::new(uf4.extract(0) * a,
                          uf4.extract(1) * a,
                          uf4.extract(2) * a,
                          uf4.extract(3) * a))
    }
}

impl Mul<f64> for Point2 {
    type Output = Point2;

    fn mul(self, a: f64) -> Point2 {
        let uf2 = self.0;
        Point2(f64x2::new(uf2.extract(0) * a, uf2.extract(1) * a))
    }
}


impl Mul<f64> for Point3 {
    type Output = Point3;

    fn mul(self, a: f64) -> Point3 {
        let uf2 = self.0;
        Point3(f64x3::new(uf2.extract(0) * a, uf2.extract(1) * a, uf2.extract(2) * a))
    }
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
impl Add for Point2 {
    type Output = Point2;

    fn add(self, p2: Point2) -> Point2 {
        Point2(self.0 + p2.0)
    }
}
impl AddAssign for Point2 {
    fn add_assign(&mut self, p2: Point2) {
        self.0 = self.0 - p2.0;
    }
}

impl Sub for Point2 {
    type Output = Point2;

    fn sub(self, p2: Point2) -> Point2 {
        Point2(self.0 - p2.0)
    }
}

impl SubAssign for Point2 {
    fn sub_assign(&mut self, p2: Point2) {
        self.0 = self.0 - p2.0;
    }
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
impl Add for Point4 {
    type Output = Point4;

    fn add(self, p2: Point4) -> Point4 {
        Point4(self.0 + p2.0)
    }
}

impl AddAssign for Point4 {
    fn add_assign(&mut self, p2: Point4) {
        self.0 = self.0 + p2.0;
    }
}

impl Sub for Point4 {
    type Output = Point4;

    fn sub(self, p2: Point4) -> Point4 {
        Point4(self.0 - p2.0)
    }
}

impl SubAssign for Point4 {
    fn sub_assign(&mut self, p2: Point4) {
        self.0 = self.0 - p2.0;
    }
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
impl Add for Point3 {
    type Output = Point3;

    fn add(self, p2: Point3) -> Point3 {
        Point3(self.0 + p2.0)
    }
}

impl AddAssign for Point3 {
    fn add_assign(&mut self, p2: Point3) {
        self.0 = self.0 + p2.0;
    }
}

impl Sub for Point3 {
    type Output = Point3;

    fn sub(self, p2: Point3) -> Point3 {
        Point3(self.0 - p2.0)
    }
}

impl SubAssign for Point3 {
    fn sub_assign(&mut self, p2: Point3) {
        self.0 = self.0 - p2.0;
    }
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
impl Point2 {
    pub fn new(x: f64, y: f64) -> Point2 {
        Point2(f64x2::new(x, y))
    }
    pub fn splat(x: f64) -> Point2 {
        Point2(f64x2::splat(x))
    }
    pub fn load(x: &[f64], i: usize) -> Point2 {
        Point2(f64x2::load(x, i))
    }
    pub fn extract(&self, i: u32) -> f64 {
        self.0.extract(i)
    }
    pub fn replace(&mut self, i: u32, x: f64) -> Point2 {
        Point2(self.0.replace(i, x))
    }
}

impl Point4 {
    pub fn new(x: f64, y: f64, z: f64, w: f64) -> Point4 {
        Point4(f64x4::new(x, y, z, w))
    }
    pub fn splat(x: f64) -> Point4 {
        Point4(f64x4::splat(x))
    }
    pub fn load(x: &[f64], i: usize) -> Point4 {
        Point4(f64x4::load(x, i))
    }
    pub fn extract(&self, i: u32) -> f64 {
        self.0.extract(i)
    }
    pub fn replace(&mut self, i: u32, x: f64) -> Point4 {
        Point4(self.0.replace(i, x))
    }
}

impl Point3 {
    pub fn new(x: f64, y: f64, z: f64) -> Point3 {
        Point3(f64x3::new(x, y, z))
    }
    pub fn splat(x: f64) -> Point3 {
        Point3(f64x3::splat(x))
    }
    pub fn load(x: &[f64], i: usize) -> Point3 {
        Point3(f64x3::load(x, i))
    }
    pub fn extract(&self, i: u32) -> f64 {
        self.0.extract(i)
    }
    pub fn replace(&mut self, i: u32, x: f64) -> Point3 {
        Point3(self.0.replace(i, x))
    }
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pub fn lerp<T: Mul<f64, Output = T> + Copy + Add<T, Output = T>>(l: f64, p: T, q: T) -> T {
    let m = 1.0 - l;
    p * m + q * l
}


pub fn dlerp<T: Mul<f64, Output = T> + Copy + Clone + SubAssign<T>>(l: f64, p: T, q: T) -> T {
    let mut r = q.clone();
    r -= p;
    r * l
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
impl Default for Point2 {
    fn default() -> Point2 {
        Point2::splat(0.0)
    }
}

impl Default for Point4 {
    fn default() -> Point4 {
        Point4::splat(0.0)
    }
}

impl Default for Point3 {
    fn default() -> Point3 {
        Point3::splat(0.0)
    }
}

pub type Point1 = f64;
// ~~~~~~~~~~~~~~~~
pub trait HiType
{
    type NextT;
}

pub trait LoType
{
    type PrevT;
}

impl HiType for Point1 {
    type NextT = Point1;
}

impl HiType for Point2 {
    type NextT = Point3;
}

impl HiType for Point3 {
    type NextT = Point4;
}
// ~~~~~~~~~~~~~~~~
impl LoType for Point2 {
    type PrevT = Point1;
}

impl LoType for Point3 {
    type PrevT = Point2;
}

impl LoType for Point4 {
    type PrevT = Point3;
}


// impl Zero for Point2
// {
//     fn zero() -> Point2 { Point2::splat(0.0) }
// }
// impl Zero for Point4
// {
//     fn zero() -> Point4 { Point4::splat(0.0) }
// }

#[test]
pub fn it_works() {
    let p = Point2::new(1.0, 2.0);
    let q: Point2 = Default::default();
    let r = lerp(0.2, p, q);
    println!("{:?}, {:?}", r, dlerp(0.2, p, q));
}
