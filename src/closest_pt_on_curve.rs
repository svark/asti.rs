// approximate distance function as a bspline and find its minima
use vectorspace::PointT;
use nalgebra::{Dot,Norm};
use spline_approx::cubic_approx1d;
use curve::{FiniteCurve,CurvePoint};
use util;
fn closest_pt_on_curve<P> (p:&P, spl: &Bspline<P>) 
 where P: PointT, <P as PointAsVector>: Dot<f64>
{
    let f = |u| -> f64 { let vec = spl.eval(u) - *p; vec.dot(&vec) };
    let to_merge = Vec::with_capacity(spl.knots().len());
    for ts in spl.knots().as_slice().windows(2) {
        if( Param(ts[0]) == Param(ts[1])) continue;
        to_merge.push(ts[0] + ts[1]);
    }
    let new_ts = util::merge(spl.knots(), to_merge)
    let quasi_interp = cubic_approx1d(f, new_ts);
    let t = rootfinder::find_next_root(quasi_interp, spl.start_param());
    // this is broken
    return CurvePoint::new(spl.eval(t), t)
}