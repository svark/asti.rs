pub use nalgebra::{Vector1, Vector2, Vector3, Vector4};
pub use nalgebra::{Point1, Point2, Point3, Point4, center, PointAsVector, Norm, Dot, Absolute,
                   zero};
pub type Pt1 = Point1<f64>;
pub type Pt2 = Point2<f64>;
pub type Pt3 = Point3<f64>;
pub type Pt4 = Point4<f64>;
pub type Vec1 = Vector1<f64>;
pub type Vec2 = Vector2<f64>;
pub type Vec3 = Vector3<f64>;
pub type Vec4 = Vector4<f64>;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use vectorspace::{PointT, VectorT, Ops};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

impl Ops for Pt1 {}
impl Ops for Pt2 {}
impl Ops for Pt3 {}
impl Ops for Pt4 {}

impl PointT for Pt1 {
    type L = Pt1;
    type H = Pt2;

    fn ldim(&self) -> Pt1 {
        panic!("hit the lower limit");
    }

    fn hdim(&self, pad: f64) -> Pt2 {
        Pt2::new(self[0], pad)
    }

    fn lerp(&self, l: f64, q: Pt1) -> Pt1 {
        let p = *self;
        p + (q - p) * l
    }
}

impl PointT for Pt2 {
    type L = Pt1;
    type H = Pt3;

    fn ldim(&self) -> Pt1 {
        Pt1::new(self[0])
    }

    fn hdim(&self, pad: f64) -> Pt3 {
        Pt3::new(self[0], self[1], pad)
    }

    fn lerp(&self, l: f64, q: Pt2) -> Pt2 {
        let p = *self;
        p + (q - p) * l
    }
}

impl PointT for Pt3 {
    type L = Pt2;
    type H = Pt4;

    fn ldim(&self) -> Pt2 {
        Pt2::new(self[0], self[1])
    }

    fn hdim(&self, pad: f64) -> Pt4 {
        Pt4::new(self[0], self[1], self[2], pad)
    }
    fn lerp(&self, l: f64, q: Pt3) -> Pt3 {
        let p = *self;
        p + (q - p) * l
    }
}


impl PointT for Pt4 {
    type L = Pt3;
    type H = Pt4;

    fn ldim(&self) -> Pt3 {
        Pt3::new(self[0], self[1], self[2])
    }

    fn hdim(&self, _: f64) -> Pt4 {
        panic!("hit the upper limit");
    }

    fn lerp(&self, l: f64, q: Pt4) -> Pt4 {
        let p = *self;
        p + (q - p) * l
    }
}

pub fn is_coplanar(p1: &Pt3, p2: &Pt3, p3: &Pt3, p4: &Pt3, tol: f64) -> bool {
    let pa = center(&center(p1, p2), &center(p3, p4));
    let v1 = *p1 - pa;
    let v2 = *p2 - pa;
    let v3 = *p3 - pa;
    let v4 = *p4 - pa;
    let (x1, y1, z1) = (v1.extract(0), v1.extract(1), v1.extract(2));
    let (x2, y2, z2) = (v2.extract(0), v2.extract(1), v2.extract(2));
    let (x3, y3, z3) = (v3.extract(0), v3.extract(1), v3.extract(2));
    let (x4, y4, z4) = (v4.extract(0), v4.extract(1), v4.extract(2));

    let det = y4 * (z3 * (x2 - x1) + z2 * (x1 - x3)) + y3 * (z4 * (x1 - x2) - x1 * z2) +
              z4 * (x3 * y2 - x1 * y2) + x1 * y2 * z3 +
              z1 * (y4 * (x3 - x2) + x2 * y3 - x3 * y2 + x4 * (y2 - y3)) +
              y1 * (z4 * (x2 - x3) - x2 * z3) + x3 * y1 * z2 +
              x4 * (-y1 * z2 + y1 * z3 - y2 * z3 + y3 * z2);

    det.abs() < tol
}

pub fn is_collinear<P: PointT>(p1: &P, p2: &P, p3: &P, tol: f64) -> bool
    where <P as PointAsVector>::Vector: Dot<f64> + Norm<NormType = f64>
{
    let p12 = *p1 - *p2;
    let p23 = *p2 - *p3;
    p12.dot(&p23).abs() < tol * p12.norm() * p23.norm()
}

impl Ops for Vec1 {}
impl Ops for Vec2 {}
impl Ops for Vec3 {}
impl Ops for Vec4 {}
// impl Lerp for Vec1{}
// impl Lerp for Vec2{}
// impl Lerp for Vec3{}
// impl Lerp for Vec4{}

impl VectorT for Vec1 {
    type L = Vec1;
    type H = Vec2;
    type P = Pt1;

    fn ldim(&self) -> Vec1 {
        panic!("hit the lower limit");
    }

    fn hdim(&self, pad: f64) -> Vec2 {
        Vec2::new(self[0], pad)
    }

    fn to_pt(&self) -> Self::P {
        self.to_point()
    }
}

impl VectorT for Vec2 {
    type L = Vec1;
    type H = Vec3;
    type P = Pt2;

    fn ldim(&self) -> Vec1 {
        Vec1::new(self[0])
    }

    fn hdim(&self, pad: f64) -> Vec3 {
        Vec3::new(self[0], self[1], pad)
    }

    fn to_pt(&self) -> Self::P {
        self.to_point()
    }
}

impl VectorT for Vec3 {
    type L = Vec2;
    type H = Vec4;
    type P = Pt3;
    fn ldim(&self) -> Vec2 {
        Vec2::new(self[0], self[1])
    }

    fn hdim(&self, pad: f64) -> Vec4 {
        Vec4::new(self[0], self[1], self[2], pad)
    }

    fn to_pt(&self) -> Self::P {
        self.to_point()
    }
}


impl VectorT for Vec4 {
    type L = Vec3;
    type H = Vec4;
    type P = Pt4;
    fn hdim(&self, _: f64) -> Vec4 {
        panic!("hit the upper limit");
    }

    fn ldim(&self) -> Vec3 {
        Vec3::new(self[0], self[1], self[2])
    }

    fn to_pt(&self) -> Self::P {
        self.to_point()
    }
}


#[test]
pub fn it_works() {
    let p = Pt2::new(1.0, 2.0);
    let q: Pt2 = Pt2::splat(0.0);
    let r = p.lerp(0.2, q);
    assert!(((r - p).norm() / (q - r).norm() - 0.25) < 1e-8);
}
