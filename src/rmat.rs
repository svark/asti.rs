use std::f64::INFINITY;
use tol::PARAMRES;
use util::upper_bound;
use tol::Tol; // for small_param
use std::ops::MulAssign;

struct RMatTau<'a, 'b> {
    knots: &'a [f64],
    taus: &'b [f64],
    offset: usize,
}

struct RMatExplicitDer<'a, 'b>(RMatTau<'a, 'b>);

trait Entries
{
    fn get_diag_from_ndiag(v: f64) -> f64;
    fn get_ndiag(&self, sz: usize, d: usize) -> f64;
    fn get_knots(&self) -> &[f64];
}

#[inline]
pub fn sdiv(x: f64, y: f64) -> f64 {
    if y.small_param() {
        0.0f64
    } else {
        x / y
    }
}

impl<'a, 'b> Entries for RMatTau<'a, 'b> {
    #[inline]
    fn get_ndiag(&self, sz: usize, d: usize) -> f64 {
        let t = self.knots;
        let idx = sz + self.offset;
        assert!(idx >= sz);
        let denom = t[idx] - t[idx - d];
        sdiv(self.taus[d] - t[idx - d], denom)
    }

    #[inline]
    fn get_diag_from_ndiag(v: f64) -> f64 {
        1.0 - v
    }

    #[inline]
    fn get_knots(&self) -> &[f64] {
        self.knots
    }
}

impl<'a, 'b> Entries for RMatExplicitDer<'a, 'b> {
    #[inline]
    fn get_ndiag(&self, sz: usize, d: usize) -> f64 {
        let t = self.0.knots;
        let idx = sz + self.0.offset;
        assert!(idx >= sz);
        let denom = t[idx] - t[idx - d];
        if denom.abs() < PARAMRES {
            INFINITY * denom.signum()
        } else {
            d as f64 / denom
        }
    }

    #[inline]
    fn get_diag_from_ndiag(v: f64) -> f64 {
        -v
    }

    #[inline]
    fn get_knots(&self) -> &[f64] {
        self.0.knots
    }
}

struct RMatExplicitResult {
    basis: Vec<f64>,
}

impl RMatExplicitResult {
    fn with_capacity(c: usize) -> RMatExplicitResult {
        let mut b = Vec::with_capacity(c);
        b.push(1.0);
        RMatExplicitResult { basis: b }
    }

    fn drain(self) -> Vec<f64> {
        self.basis
    }
}

impl<'a, T: Entries> MulAssign<&'a T> for RMatExplicitResult {
    fn mul_assign(&mut self, mat: &'a T) {
        let basis = &mut self.basis;
        let num_cols = basis.len();
        let mut c = mat.get_ndiag(num_cols, num_cols);
        let e = basis[num_cols - 1] * c;
        basis.resize(num_cols + 1, e);
        for k in (1..num_cols).rev() {
            let b = T::get_diag_from_ndiag(c);
            c = mat.get_ndiag(k, num_cols);
            basis[k] = b * basis[k] + c * basis[k - 1];
        }
        basis[0] *= T::get_diag_from_ndiag(c);
    }
}

pub fn locate_nu(u: f64, d: usize, t: &[f64]) -> usize {
    let ncpts = t.len() - d - 1;
    let mut u = u;
    if u <= t[d] {
        u = t[d] + PARAMRES / 2.0;
    }
    if u >= t[ncpts] {
        u = t[ncpts] - PARAMRES / 2.0;
    }
    let idx = upper_bound(t, u) - 1;
    assert!(idx >= d && idx < ncpts);
    idx
}


pub fn eval(knots: &[f64], us_nu: (&[f64], Option<usize>), deg: u32, der_order: u32) -> Vec<f64> {

    let (us, nu) = us_nu;
    assert!(us.len() > 0);
    let mu = {
        if let Some(mu) = nu {
            mu
        } else {
            locate_nu(us[0], deg as usize, knots)
        }
    };
    assert!(mu < knots.len());

    let rmat_tau = RMatTau {
        knots: knots,
        taus: us,
        offset: mu,
    };

    let mut res = RMatExplicitResult::with_capacity(deg as usize + 1);

    for _ in 1..(deg - der_order + 1) {
        res *= &rmat_tau;
    }

    let rmat_der = RMatExplicitDer(rmat_tau);
    for j in (deg - der_order + 1)..(deg + 1) {
        assert_eq!(res.basis.len(), j as usize);
        res *= &rmat_der;
    }
    res.drain()
}

pub fn basis(knots: &[f64], u_nu: (f64, Option<usize>), deg: u32) -> Vec<f64> {
    eval(knots,
         (vec![u_nu.0;deg as usize + 1].as_slice(), u_nu.1),
         deg,
         0)
}
