#[cfg(target_feature = "sse2")]
use simd::x86::sse2::f64x2;

// #[cfg(target_feature = "avx")]
// use simd::x86::avx::f64x4;

// #[cfg(target_feature = "neon")]
// use simd::aarch64::neon::f64x2;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#[cfg(all(not(target_feature = "sse2"), not(target_feature = "neon")))]
#[derive(Debug, Copy, Clone)]
#[allow(non_camel_case_types)]
pub struct f64x2(f64, f64);

//#[cfg(not(target_feature = "avx"))]
#[allow(non_camel_case_types)]
#[derive(Debug, Copy, Clone)]
pub struct f64x4(f64, f64, f64, f64);

#[allow(non_camel_case_types)]
#[derive(Debug, Copy, Clone)]
pub struct f64x3(f64, f64, f64);

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#[cfg(not(target_feature = "sse2"))]
impl f64x2 {
    pub fn new(x: f64, y: f64) -> f64x2 {
        f64x2(x, y)
    }

    fn load(x: &[f64], i: usize) -> f64x2 {
        f64x2(x[i], x[i + 1])
    }

    fn splat(x: f64) -> f64x2 {
        f64x2(x, x)
    }

    fn extract(&self, i: u32) -> f64 {
        match i {
            0 => self.0,
            1 => self.1,
            _ => panic!("bad index for extract in f64x2"),
        }
    }

    fn replace(&mut self, i: u32, x: f64) -> f64x2 {
        match i {
            0 => f64x2(x, self.1),
            1 => f64x2(self.0, x),
            _ => panic!("bad index for replace in f64x2"),
        }
    }
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//#[cfg(not(target_feature = "avx"))]
impl f64x4 {
    pub fn new(x: f64, y: f64, z: f64, w: f64) -> f64x4 {
        f64x4(x, y, z, w)
    }

    fn load(x: &[f64], i: usize) -> f64x4 {
        f64x4(x[i], x[i + 1], x[i + 2], x[i + 3])
    }

    fn splat(x: f64) -> f64x4 {
        f64x4(x, x, x, x)
    }

    fn extract(&self, i: u32) -> f64 {
        match i {
            0 => self.0,
            1 => self.1,
            2 => self.2,
            3 => self.3,
            _ => panic!("bad index for extract in f64x4"),
        }
    }

    fn replace(&mut self, i: u32, x: f64) -> f64x4 {
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

    fn load(x: &[f64], i: usize) -> f64x3 {
        f64x3(x[i], x[i + 1], x[i + 2])
    }

    fn splat(x: f64) -> f64x3 {
        f64x3(x, x, x)
    }

    fn extract(&self, i: u32) -> f64 {
        match i {
            0 => self.0,
            1 => self.1,
            2 => self.2,
            _ => panic!("bad index for extract in f64x3"),
        }
    }

    fn replace(&mut self, i: u32, x: f64) -> f64x3 {
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

// #[cfg(not(target_feature = "avx"))]
impl Add for f64x4 {
    type Output = f64x4;

    fn add(self, p2: f64x4) -> f64x4 {
        f64x4(self.0 + p2.0, self.1 + p2.1, self.2 + p2.2, self.3 + p2.3)
    }
}

// #[cfg(not(target_feature = "avx"))]
impl Sub for f64x4 {
    type Output = f64x4;

    fn sub(self, p2: f64x4) -> f64x4 {
        f64x4(self.0 - p2.0, self.1 - p2.1, self.2 - p2.2, self.3 - p2.3)
    }
}

// #[cfg(not(target_feature = "avx"))]
impl Mul<f64x4> for f64x4 {
    type Output = f64x4;

    fn mul(self, p2: f64x4) -> f64x4 {
        f64x4(self.0 * p2.0, self.1 * p2.1, self.2 * p2.2, self.3 * p2.3)
    }
}

// #[cfg(not(target_feature = "avx"))]
impl Mul<f64> for f64x4 {
    type Output = f64x4;

    fn mul(self, p2: f64) -> f64x4 {
        f64x4(self.0 * p2, self.1 * p2, self.2 * p2, self.3 * p2)
    }
}



// #[cfg(not(target_feature = "avx"))]
impl AddAssign for f64x4 {
    fn add_assign(&mut self, p2: f64x4) {
        *self = *self + p2;
    }
}


// #[cfg(not(target_feature = "avx"))]
impl SubAssign for f64x4 {
    fn sub_assign(&mut self, p2: f64x4) {
        *self = *self - p2;
    }
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#[cfg(not(target_feature = "sse2"))]
impl Add for f64x2 {
    type Output = f64x2;

    fn add(self, p2: f64x2) -> f64x2 {
        f64x2(self.0 + p2.0, self.1 + p2.1)
    }
}

#[cfg(not(target_feature = "sse2"))]
impl Sub for f64x2 {
    type Output = f64x2;

    fn sub(self, p2: f64x2) -> f64x2 {
        f64x2(self.0 - p2.0, self.1 - p2.1)
    }
}

#[cfg(not(target_feature = "sse2"))]
impl Mul<f64x2> for f64x2 {
    type Output = f64x2;

    fn mul(self, p2: f64x2) -> f64x2 {
        f64x2(self.0 * p2.0, self.1 * p2.1)
    }
}

#[cfg(not(target_feature = "sse2"))]
impl Mul<f64> for f64x2 {
    type Output = f64x2;

    fn mul(self, p2: f64) -> f64x2 {
        f64x2(self.0 * p2, self.1 * p2)
    }
}



#[cfg(not(target_feature = "sse2"))]
impl AddAssign for f64x2 {
    fn add_assign(&mut self, p2: f64x2) {
        *self = *self + p2;
    }
}


#[cfg(not(target_feature = "sse2"))]
impl SubAssign for f64x2 {
    fn sub_assign(&mut self, p2: f64x2) {
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

impl Mul<f64x3> for f64x3 {
    type Output = f64x3;

    fn mul(self, p2: f64x3) -> f64x3 {
        f64x3(self.0 * p2.0, self.1 * p2.1, self.2 * p2.2)
    }
}

impl Mul<f64> for f64x3 {
    type Output = f64x3;

    fn mul(self, p2: f64) -> f64x3 {
        f64x3(self.0 * p2, self.1 * p2, self.2 * p2)
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
}

impl Point4 {
    pub fn new(x: f64, y: f64, z: f64, w: f64) -> Point4 {
        Point4(f64x4::new(x, y, z, w))
    }
}

impl Point3 {
    pub fn new(x: f64, y: f64, z: f64) -> Point3 {
        Point3(f64x3::new(x, y, z))
    }
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use vectorspace::VectorSpace;

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

impl VectorSpace for Point1 {
    fn dim() -> u32 {
        1
    }

    fn load(x: &[f64], i: usize) -> Self {
        assert_eq!(i, 0);
        x[0]
    }

    fn splat(x: f64) -> Point1 {
        x
    }

    fn extract(&self, i: u32) -> f64 {
        assert_eq!(i, 0);
        *self
    }

    fn replace(&mut self, i: u32, x: f64) -> Point1 {
        assert_eq!(i, 0);
        x
    }
    type L = Point1;
    type H = Point2;

    fn ldim(&self) -> Point1 {
        panic!("hit the lower limit");
    }

    fn hdim(&self, pad: f64) -> Point2 {
        Point2::new(*self, pad)
    }

    fn dot(&self,v:&Self) -> f64 {
        self * v
    }

    fn len(&self) -> f64 { self.abs() }
}

impl VectorSpace for Point2 {

    fn dim() -> u32 {
        2
    }

    fn splat(x: f64) -> Point2 {
        Point2(f64x2::splat(x))
    }

    fn load(x: &[f64], i: usize) -> Point2 {
        Point2(f64x2::load(x, i))
    }

    fn extract(&self, i: u32) -> f64 {
        self.0.extract(i)
    }

    fn replace(&mut self, i: u32, x: f64) -> Point2 {
        Point2(self.0.replace(i, x))
    }

    type L = Point1;
    type H = Point3;

    fn ldim(&self) -> Point1 {
        self.extract(0)
    }
    fn hdim(&self, pad: f64) -> Point3 {
        Point3::new(self.extract(0), self.extract(1), pad)
    }

    fn dot(&self,v:&Self) -> f64 {
        let u = self.0 * v.0;
        u.extract(0) + u.extract(1)
    }

}
impl VectorSpace for Point3 {
    fn dim() -> u32 {
        3
    }
    fn splat(x: f64) -> Point3 {
        Point3(f64x3::splat(x))
    }
    fn load(x: &[f64], i: usize) -> Point3 {
        Point3(f64x3::load(x, i))
    }
    fn extract(&self, i: u32) -> f64 {
        self.0.extract(i)
    }
    fn replace(&mut self, i: u32, x: f64) -> Point3 {
        Point3(self.0.replace(i, x))
    }

    type L = Point2;
    type H = Point4;

    fn ldim(&self) -> Point2 {
        Point2::new(self.extract(0), self.extract(1))
    }
    fn hdim(&self, pad: f64) -> Point4 {
        Point4::new(self.extract(0), self.extract(1), self.extract(2), pad)
    }

    fn dot(&self,v:&Self) -> f64 {
        let u = self.0 * v.0;
        u.extract(0) + u.extract(1) + u.extract(2)
    }
}

pub fn is_coplanar(p1:&Point3, p2:&Point3, p3:&Point3, p4:&Point3, tol: f64) -> bool
{
    let pa = (*p1+*p2+*p3+*p4)*0.25;
    let p1 = *p1 - pa;
    let p2 = *p2 - pa;
    let p3 = *p3 - pa;
    let p4 = *p4 - pa;
    let (x1,y1,z1) = (p1.extract(0),p1.extract(1), p1.extract(2) );
    let (x2,y2,z2) = (p2.extract(0),p2.extract(1), p2.extract(2) );
    let (x3,y3,z3) = (p3.extract(0),p3.extract(1), p3.extract(2) );
    let (x4,y4,z4) = (p4.extract(0),p4.extract(1), p4.extract(2) );

    let det = y4*(z3*(x2-x1)+z2*(x1-x3))+y3*(z4*(x1-x2)-x1*z2)+z4*(x3*y2-x1*y2)+x1*y2*z3+z1*(y4*(x3-x2)+x2*y3-x3*y2+x4*(y2-y3))+y1*(z4*(x2-x3)-x2*z3)+x3*y1*z2+x4*(-y1*z2+y1*z3-y2*z3+y3*z2);

    det.abs() < tol
}

pub fn  is_collinear<P:VectorSpace>(p1: &P, p2: &P, p3: &P, tol: f64) -> bool
{
    let p12 = *p1 - *p2;
    let p23 = *p2 - *p3;
    p12.dot(&p23).abs() - p12.len() * p23.len() < tol*tol
}

impl VectorSpace for Point4 {
    fn dim() -> u32 {
        4
    }
    fn splat(x: f64) -> Point4 {
        Point4(f64x4::splat(x))
    }
    fn load(x: &[f64], i: usize) -> Point4 {
        Point4(f64x4::load(x, i))
    }
    fn extract(&self, i: u32) -> f64 {
        self.0.extract(i)
    }
    fn replace(&mut self, i: u32, x: f64) -> Point4 {
        Point4(self.0.replace(i, x))
    }
    type L = Point3;
    type H = Point4;

    fn ldim(&self) -> Point3 {
        Point3::new(self.extract(0), self.extract(1), self.extract(2))
    }

    fn hdim(&self, _: f64) -> Point4 {
        panic!("hit the upper limit")
    }

    fn dot(&self,v:&Self) -> f64 {
        let u = self.0 * v.0;
        u.extract(0) + u.extract(1) + u.extract(2) + u.extract(3)
    }
}


#[test]
pub fn it_works() {
    let p = Point2::new(1.0, 2.0);
    let q: Point2 = Default::default();
    let r = p.lerp(0.2, q);
    assert!(((r-p).len()/(q-r).len() - 0.25) < 1e-8);
}
