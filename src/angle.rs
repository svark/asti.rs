use point::{Pt3, Vec2, Vec3};
use nalgebra::{RotationTo, Cross, Norm, FloatVector};
use tol::{Tol, RESABS};
use vectorspace::Ops;

pub trait Angle : FloatVector<f64> + RotationTo
{
    fn angle(&self, other: &Self) -> <Self as RotationTo>::AngleType {
        self.angle_to(other)
    }

    fn try_angle(&self, other: &Self) -> Option<<Self as RotationTo>::AngleType> {
        if self.norm().small() || other.norm().small() {
            return None;
        }
        return Some(self.angle(other));
    }

    fn angle_between_unitvecs(&self, other: &Self) -> f64 {
        self.dot(&other).acos()
    }
}


impl Angle for Vec3 {}
impl Angle for Vec2 {}

pub fn plane_normal(p: [Pt3; 3]) -> Option<Vec3> {
    let ref v10 = p[1] - p[0];
    let ref v20 = p[2] - p[0];
    v10.cross(&v20).try_normalize(RESABS)
}

pub fn perp_in_plane(v1: Vec3, p: [Pt3; 3]) -> Option<Vec3> {
    if let Some(v) = plane_normal(p) {
        Some(v1.cross(&v))
    } else {
        None
    }
}

pub fn perp_in_2dplane(v1: Vec2) -> Vec2 {
    Vec2::new(-v1.extract(1), v1.extract(0))
}
