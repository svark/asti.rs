pub const PARAMRES: f64 = 1e-10;
/// parameter resolution
pub const RESABS: f64 = 1e-8;
/// model resolution
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

#[derive(Debug,Clone,Copy)]
pub struct Param(pub f64);

pub fn wrap_param(v: f64) -> Param {
    return Param(v);
}
pub fn unwrap_param(s: Param) -> f64 {
    let Param(v) = s;
    v
}
pub fn unwrap_param_ref(s: &Param) -> f64 {
    let &Param(v) = s;
    v
}

impl PartialEq for Param {
    fn eq(&self, other: &Self) -> bool {
        (self.0 - other.0).small_param()
    }
}

impl Eq for Param {}

impl PartialOrd for Param {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.eq(other) {
            Some(Ordering::Equal)
        } else {
            self.0.partial_cmp(&other.0)
        }
    }
}
// workaround for lack of total order in floating point values (courtesy NANs)..
impl Ord for Param {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

#[derive(Debug,Clone,Copy)]
pub struct OSCoord(pub f64);

pub fn wrap_oscoord(v: f64) -> OSCoord {
    return OSCoord(v);
}

pub fn unwrap_oscoord(s: OSCoord) -> f64 {
    let OSCoord(v) = s;
    v
}
pub fn unwrap_oscoord_ref(s: &OSCoord) -> f64 {
    let &OSCoord(v) = s;
    v
}

impl PartialEq for OSCoord {
    fn eq(&self, other: &Self) -> bool {
        (self.0 - other.0).small()
    }
}

impl Eq for OSCoord {}

impl PartialOrd for OSCoord {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.eq(other) {
            Some(Ordering::Equal)
        } else {
            self.0.partial_cmp(&other.0)
        }
    }
}
// workaround for lack of total order in floating point values (courtesy NANs)..
impl Ord for OSCoord {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}
