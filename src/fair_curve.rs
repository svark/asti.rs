use tol::{PARAMRES, Param};
use bspline::Bspline;
use errcodes::GeomErrorCode;
use vectorspace::{PointT, to_pt};
use curve::Curve;
use splinedata::{KnotManip, SplineData};
use itertools::Itertools;
use nalgebra::Norm;

// Farin, Curves and Surfaces pg 422, Shape
// Idea is to move a control point where the largest jump is noticed in third derivatives.
// control points corresponding to these points then moved to a more 'favourable' location
// where the offending knot may be deemed removable
pub fn fair_curve<P>(crv: Bspline<P>, tol: f64) -> Result<Bspline<P>, GeomErrorCode>
    where P: PointT
{
    if crv.degree() < 3 {
        return Err(GeomErrorCode::GeomDegreeLow);
    }
    let ts = crv.knots();
    let tlen = ts.len();
    let d = crv.degree() as usize;
    let uts: Vec<f64> = uniq_ts!(ts[d..tlen - d].iter()).collect();

    let mut kdash_variation = Vec::with_capacity(uts.len());
    for u in uts.iter() {
        let kdash_minus = crv.eval_derivative(u - PARAMRES / 2.0, 3);
        let kdash_plus = crv.eval_derivative(u + PARAMRES / 2.0, 3);
        kdash_variation.push((kdash_plus - kdash_minus).norm());
    }
    debug_assert!(kdash_variation.len() > 0);
    let (k, _) = kdash_variation.iter()
                                .enumerate()
                                .max_by_key(|&(_, a)| Param(*a))
                                .expect("empty");

    let u = uts[k];
    let nu = crv.locate_nu(u);
    let cpts = crv.control_points();
    let t = ts;
    let mut newcpts = cpts.clone();

    let mut l = cpts[nu - 1].to_vector() * (t[nu + 1] - t[nu - 3]) -
                cpts[nu - 2].to_vector() * (t[nu + 1] - t[nu]);

    l *= 1.0 / (t[nu] - t[nu - 3]);

    let mut r = (cpts[nu + 1] * (t[nu + 3] - t[nu - 1])).to_vector() -
                (cpts[nu + 2] * (t[nu] - t[nu - 1])).to_vector();

    r *= 1.0 / (t[nu] - t[nu - 3]);

    newcpts[nu] = to_pt(l * (t[nu + 2] - t[nu]) + r * (t[nu] - t[nu - 2]));
    newcpts[nu] *= 1.0 / (t[nu + 2] - t[nu - 2]);

    let den = (newcpts[nu] - cpts[nu]).norm();
    if den > tol {
        newcpts[nu] = cpts[nu] + (newcpts[nu] - cpts[nu]) * tol / den.sqrt();
    } else {
        newcpts[nu] = cpts[nu];
    }

    let newts = ts.clone();
    return Ok(Bspline::new(newcpts, newts));
}
