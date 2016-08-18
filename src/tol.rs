pub const PARAMRES: f64 = 1e-10;
pub const RESABS: f64 = 1e-8;
pub const SQRESABS: f64 = 1e-16;
use std::cmp::{Ord, Ordering};

pub trait Tol
{
    fn eqres(self, w: f64) -> bool;
    fn eqparam(self, w: f64) -> bool;
    fn small(self) -> bool;
    fn small_param(self) -> bool;
}

impl Tol for f64 {
    #[inline]
    fn eqres(self, w: f64) -> bool {
        return (self - w).abs() < RESABS;
    }

    #[inline]
    fn eqparam(self, w: f64) -> bool {
        return (self - w).abs() < PARAMRES;
    }

    #[inline]
    fn small(self) -> bool {
        return self.abs() < RESABS;
    }

    #[inline]
    fn small_param(self) -> bool {
        return self.abs() < PARAMRES;
    }
}

pub struct Param(pub f64);

pub fn wrap_param(v: f64) -> Param { return Param(v); }
pub fn unwrap_param(s: Param) -> f64 { let Param(v) = s; v }
pub fn unwrap_param_ref(s: &Param) -> f64 { let &Param(v) = s; v }

impl PartialEq for Param
{
    fn eq( &self, other: &Self) -> bool
    {
        (self.0 - other.0).small_param()
    }
}

impl Eq for Param {}

impl PartialOrd for Param
{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.eq(other) {
            Some(Ordering::Equal)
        }else {
            self.0.partial_cmp(&other.0)
        }
    }
}
// workaround for lack of total order in floating point values (courtesy NANs)..
impl Ord for Param
{
    fn cmp(&self, other: &Self) -> Ordering {
         self.partial_cmp(other).unwrap()
    }

}
