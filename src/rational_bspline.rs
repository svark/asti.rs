use point::{HiType};
use bspline::Bspline;

pub struct RationalBspline<Point:HiType>  {
    pub spl: Bspline<Point::NextT>,
}
