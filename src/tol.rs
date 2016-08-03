pub const PARAMRES: f64 = 1e-10;
pub const RESABS: f64 = 1e-8;
pub const SQRESABS: f64 = 1e-16;

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
