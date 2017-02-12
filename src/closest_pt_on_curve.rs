// approximate distance function as a bspline and find its minima
use vectorspace::PointT;
use nalgebra::Dot;
use spline_approx::cubic_approx1d;
use curve::{CurvePoint, Domain};
use bspline::Bspline;
use tol::{Param, Tol};
use splinedata::{SplineData, KnotManip};
use curve::Curve;
use box_compute::find_min_bound_by_insertion;
use vectorspace::Ops;

pub fn closest_pt_on_curve<P>(p: &P, spl: &Bspline<P>) -> CurvePoint<P>
    where P: PointT
{
    let f = |u: f64| {
        let vec = spl.eval(u) - *p;
        vec.dot(&vec)
    };
    let mut newts: Vec<f64> = Vec::with_capacity(2 * spl.knots().len());
    for _ in 0..4 {
        newts.push(spl.front());
    }
    for ts in spl.knots().as_slice().windows(2) {
        if Param(ts[0]) != Param(ts[1]) {
            newts.push(0.5 * (ts[0] + ts[1]));
            newts.push(ts[1]);
        }
    }
    let lastt = spl.back();
    for _ in 0..3 {
        newts.push(lastt);
    }

    let quasi_interp = cubic_approx1d(&f, newts);
    // assert!((quasi_interp.eval(-1.0).extract(0) - f(-1.0)).abs().small());
    assert!((quasi_interp.eval(0.0).extract(0) - f(0.0)).abs().small());
    assert!((quasi_interp.eval(1.0).extract(0) - f(1.0)).abs().small());

    let (t, _) = find_min_bound_by_insertion(&quasi_interp, 0);
    let fdiv = |u: f64| {
        let vec = spl.eval(u) - *p;
        let der1 = spl.eval_derivative(u, 1);
        let vecdash = der1.as_vector();
        2.0 * vec.dot(vecdash)
    };
    let fdivdiv = |u: f64| {
        let vec = spl.eval(u) - *p;
        let ref der1 = spl.eval_derivative(u, 1);
        let vecdash = der1.as_vector();
        let ref der2 = spl.eval_derivative(u, 2);
        let vecdashdash = der2.as_vector();
        2.0 * vec.dot(vecdashdash) + 2.0 * vecdash.dot(vecdash)
    };

    let t_at_min_dist = {
        let mut num_iter = 0;
        let mut u = t;
        let (s, e) = spl.param_range();
        while !fdiv(t).abs().small() && !fdivdiv(t).abs().small() && num_iter < 5 && u < e &&
              u > s {
            u = u - fdiv(u) / fdivdiv(u);
            num_iter += 1;
        }
        u
    };
    CurvePoint {
        pnt: spl.eval(t_at_min_dist),
        t: t_at_min_dist,
    }
}

#[test]
fn it_works() {
    use monomial_form::MonomialForm;
    use point::Pt1;
    let mf = MonomialForm::new(vec![Pt1::new(1.0), Pt1::new(-2.0), Pt1::new(4.)], -1.0, 1.0);
    use change_basis::to_bezier;
    let bzf = to_bezier(&mf);
    println!("at 0:{:?}", mf.eval(0.0).extract(0));
    assert!((bzf.eval(0.0).extract(0) - mf.eval(0.0).extract(0)).abs().small());
    assert!((bzf.eval(0.0).extract(0) - 1.0).abs().small());
    let cp = closest_pt_on_curve(&Pt1::new(0.0), &*bzf);
    assert!((cp.t + 0.5).small());
}