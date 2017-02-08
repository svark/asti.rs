use errcodes::GeomErrorCode;
use bspline::Bspline;
use splinedata::{SplineData, KnotManip};
use nalgebra::Norm;
use smat::{clamp_at_left, rebase_at_left, rebase_at_right};
use raise_degree::match_degrees;
use vectorspace::PointT;
use std::mem::swap;
use std::iter::once;
use rev::reverse_curve;
use curve::Curve;
use reparametrize::{reparametrize_start, reparametrize};
use split_curve::{split_open_curve, split_periodic_curve};

pub fn join_starts<T: PointT>(spl1: &Bspline<T>,
                              spl2: &Bspline<T>,
                              join_cont: usize)
                              -> Result<Bspline<T>, GeomErrorCode> {
    let mut spl1_clamped = clamp_at_left(spl1.front(), spl1);
    let mut spl2_clamped = clamp_at_left(spl2.front(), spl2);

    match_degrees(&mut spl1_clamped, &mut spl2_clamped);

    let p = spl1_clamped.degree() as usize;
    if join_cont >= p {
        return Err(GeomErrorCode::JoinContinuityTooTight);
    }
    // let the two splines have the same start parameters
    let spl1_clamped = {
        reparametrize_start(spl1_clamped, 0.0)
    };
    let spl2_clamped = {
        reparametrize_start(spl2_clamped, 0.0)
    };

    let mut a = spl1_clamped.knots()[p + 1];
    let mut b = spl2_clamped.knots()[p + 1];

    if a > b {
        swap(&mut a, &mut b)
    }
    let mut ks: Vec<f64> = vec![-a;join_cont + 1];
    ks.extend(vec![0.0;p - join_cont]);

    let enrich_knots_at_a = |spl: Bspline<T>| {
        let mut nk: i32 = (join_cont + 1) as i32;
        nk -= spl.mult(a) as i32;
        if nk > 0 {
            let as_: Vec<f64> = vec![a;nk as usize];
            spl.insert_knots(&as_)
        } else {
            spl
        }
    };
    let spl2_clamped = enrich_knots_at_a(spl2_clamped);
    let spl1_clamped = enrich_knots_at_a(spl1_clamped);

    let spl1_clamped = {
        rebase_at_left(&spl1_clamped, 0.0, ks.as_slice())
    };
    let spl2_clamped = {
        rebase_at_left(&spl2_clamped, 0.0, ks.as_slice())
    };

    let cpts1 = spl1_clamped.control_points();
    let cpts2 = spl2_clamped.control_points();

    let mut cpts: Vec<T> = Vec::with_capacity(cpts1.len() + cpts2.len());
    for &cpt in cpts1.iter().skip(join_cont + 1).rev() {
        cpts.push(cpt);
    }

    for (cpt1, cpt2) in cpts1[..join_cont + 1]
                            .into_iter()
                            .rev()
                            .zip(cpts2[..join_cont + 1].into_iter()) {
        cpts.push(cpt1.lerp(0.5, *cpt2))
    }

    for &cpt2 in cpts2[join_cont + 1..].iter() {
        cpts.push(cpt2)
    }

    // merge the two knot sequences together to get the knot sequence
    // for the join
    let s1 = spl1_clamped.knots().len();
    let spl1_clamped = {
        reverse_curve(&spl1_clamped)
    };

    let newknots: Vec<f64> = spl1_clamped.knots()[0..s1 - (join_cont + 1)]
                                 .iter()
                                 .cloned()
                                 .chain(spl2_clamped.knots()[(p + 1)..].iter().cloned())
                                 .collect();

    Ok(Bspline::new(cpts, newknots))
}

pub fn extend_curve_end_to_pt<T>(spl: &Bspline<T>, target: &T) -> Bspline<T>
    where T: PointT
{
    let clmped: Bspline<T> = clamp_at_left(spl.front(), spl);
    let s = reparametrize(clmped, 0., 1.);
    let t = s.knots();
    let n = s.control_points().len();
    let d = s.degree() as usize;
    let mut chord_len = 0.0f64;
    let mut pt = s.eval(t[d]);

    for r in 0..n - d {
        let newc = s.eval(t[d + r + 1]);
        chord_len += (newc - pt).norm();
        pt = newc;
    }

    let ld = (*target - pt).norm();
    chord_len += ld;
    let delta = ld / chord_len;

    let mut ks = vec![1.0 + delta; d + 1];
    ks[0] = 1.0;

    let exs = rebase_at_right(&s, s.back(), ks.as_slice());
    let newks: Vec<f64> = exs.knots().iter().cloned().chain(once(1.0 + delta)).collect();
    let newcpts = exs.control_points().iter().cloned().chain(once(*target)).collect();

    Bspline::new(newcpts, newks)
}
use curve::FiniteCurve;
use periodic_spline::{periodic_param, PeriodicBspline};

pub fn trim_open_curve<T: PointT>(spl: &Bspline<T>, a: f64, b: f64) -> Bspline<T> {
    assert!(b >= a);
    assert!(spl.start_param() <= a);
    assert!(spl.end_param() >= b);
    let (_, sa) = split_open_curve(&spl, a);
    let (sb, _) = split_open_curve(&sa, b);
    return sb;
}

pub fn trim_periodic_curve<T: PointT>(spl: &PeriodicBspline<T>, a: f64, b: f64) -> Bspline<T> {
    let a = periodic_param(spl.param_range(), a);
    let b = periodic_param(spl.param_range(), b);
    let sa = split_periodic_curve(&spl, a);
    let bbar = sa.start_param() + (b - spl.start_param());
    let (sb, _) = split_open_curve(&sa, bbar);
    return sb;
}


#[test]
fn it_works() {
    use point::Pt1;
    use vectorspace::Ops;
    use tol::Tol;
    let bspl = Bspline::new(vec![Pt1::new(0.0), Pt1::new(0.5), Pt1::new(0.2)],
                            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
    let bspl_ex = extend_curve_end_to_pt(&bspl, &Pt1::new(0.05));
    assert!(bspl_ex.control_points().last().unwrap().extract(0).eqres(0.05));
    let p = bspl_ex.eval(0.1);
    use closest_pt_on_curve::closest_pt_on_curve;
    let crvpt = closest_pt_on_curve(&p, &bspl);
    assert!((crvpt.pnt - p).norm() < 1e-6);

    let (sa, sb) = split_open_curve(&bspl, 0.7);
    use curve::Curve;
    use rev::reverse_curve;
    let rsa = reverse_curve(&sa);
    let joined_curve = join_starts(&rsa, &sb, 1).unwrap();
    assert!((joined_curve.eval(0.1) - bspl.eval(0.8)).norm().small());
}
