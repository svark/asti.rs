use bspline::Bspline;
use periodic_spline::{PeriodicBspline, periodic_param};
use tol::Tol;
use smat::{clamp_ends, clamp_at_right, clamp_at_left};
use splinedata::KnotManip;
use rmat::locate_nu;
use vectorspace::PointT;
use splinedata::SplineData;
use curve::Domain;

pub fn split_open_curve<P: PointT>(spl: &Bspline<P>, u: f64) -> (Bspline<P>, Bspline<P>) {
    (clamp_at_right(u, spl), clamp_at_left(u, spl))
}

pub fn rotate_base_knot<P: PointT>(pspl: &PeriodicBspline<P>, nu: usize) -> PeriodicBspline<P> {
    let p = pspl.degree() as usize;
    assert!(nu >= p);
    let cpts = pspl.control_points();
    let rng0 = pspl.param_range();
    let t = pspl.knots();
    let cpts_unwrapped = cpts[p..].to_owned();
    let sz = cpts_unwrapped.len();
    let cpts_unwrapped = cpts_unwrapped.iter()
                                       .clone()
                                       .cycle()
                                       .skip(nu - p)
                                       .take(sz)
                                       .map(|&x| x)
                                       .collect::<Vec<_>>();
    let rng = (t[nu], t[nu] + rng0.1 - rng0.0);
    let mut knots_unwrapped: Vec<f64> = t[p..t.len() - p]
                                            .into_iter()
                                            .clone()
                                            .cycle()
                                            .skip(nu - p)
                                            .take(t.len() - 2 * p - 1)
                                            .map(|&x| periodic_param(rng, x))
                                            .collect();

    let mut ps = periodic_param(rng, t[t.len() - p - 1]);
    for chks in t[p..nu + 1].windows(2) {
        ps += chks[1] - chks[0];
        knots_unwrapped.push(ps);
    }
    PeriodicBspline::new_from_unwrapped(cpts_unwrapped, knots_unwrapped, p)
}

pub fn split_periodic_curve<P: PointT>(pspl: &PeriodicBspline<P>, u: f64) -> Bspline<P> {
    let p = pspl.degree() as usize;
    if periodic_param(pspl.param_range(), u).eqparam(pspl.start_param()) {
        let cbspl = pspl.as_ref().clone();
        return clamp_ends(cbspl);
    }

    let ts = pspl.knots();
    let nu = locate_nu(u, p, ts);
    if !ts[nu].eqparam(u) {
        return split_periodic_curve(&pspl.insert_knot(u), u);
    }
    let pc = rotate_base_knot(pspl, nu);
    let cbspl = pc.into();
    return clamp_ends(cbspl);
}

#[test]
fn it_works() {
    use point::Pt2;
    use curve::Curve;
    use nalgebra::Norm;
    let pts = vec![Pt2::new(0.0, 0.0), Pt2::new(0.4, 0.3), Pt2::new(0.2, 0.8), Pt2::new(-0.2, 0.4)];
    let ks = vec![0.0, 0.3, 0.6, 0.8, 1.0];

    let pbs = PeriodicBspline::new_from_unwrapped(pts, ks, 2);
    let spbs = split_periodic_curve(&pbs, 0.2);
    assert!((spbs.eval(0.35) - pbs.eval(0.35)).norm().small());
}