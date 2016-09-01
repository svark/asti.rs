use  vectorspace::VectorSpace;
use point::{Point2, Point3};
pub trait Angle : VectorSpace
{
    fn angle(&self, other: &Self) -> Option<f64> {
        if let (Some(u1),Some(u2)) = (self.normalize(),other.normalize())
        {
            Some(u1.dot(&u2).acos())
        }else {
            None
        }
    }

    fn angle_between_unitvecs(&self, other: &Self) -> f64 {
        debug_assert!(self.angle(other).is_some());
        self.dot(&other).acos()
    }
}

impl Angle for Point3{}
impl Angle for Point2{}

pub fn cross(v1: &Point3, v2: &Point3) -> Point3
{
    let v3 = Point3::new(
        (v1.extract(1)*v2.extract(2) - v1.extract(2)*v2.extract(1)),
        (v1.extract(2)*v2.extract(0) - v1.extract(0)*v2.extract(2)),
        (v1.extract(0)*v2.extract(1) - v1.extract(1)*v2.extract(0)));
    v3
}

pub fn plane_normal(p : [Point3;3]) -> Option<Point3>
{
    let ref p10 = p[1] - p[0];
    let ref p20 = p[2] - p[0];
    cross( p10, p20 ).normalize()
}

pub fn perp_in_plane(v1: Point3, p: [Point3;3]) -> Option<Point3>
{
    if let Some(v) = plane_normal(p)
    {
        Some(cross(&v1, &v))
    }else
    {
        None
    }
}

pub fn perp_in_2dplane(v1: Point2) -> Point2
{
    Point2::new(-v1.extract(1), v1.extract(0))
}
