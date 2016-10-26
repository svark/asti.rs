use errcodes::GeomErrorCode;
use bspline::Bspline;
use splinedata::{SplineData, KnotManip};
use nalgebra::{Norm, PointAsVector};
use smat::{clamp_at_left, rebase_at_left, rebase_at_right};
use raise_degree::match_degrees;
use vectorspace::PointT;
use std::mem::swap;
use std::iter::{repeat, once};
use rev::reverse_curve;
use curve::Curve;
use reparametrize::{reparametrize_start, reparametrize};

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
    let ks: Vec<f64> = repeat(-a)
                           .take(join_cont + 1)
                           .chain(repeat(0.0).take(p - join_cont))
                           .collect();

    let nk = join_cont + 1 - spl2_clamped.mult(a);

    let as_: Vec<f64> = vec![a;nk];
    let spl2_clamped = {
        spl2_clamped.insert_knots(&as_)
    };

    let nk = join_cont + 1 - spl1_clamped.mult(a);
    let as_: Vec<f64> = vec![a;nk];
    let spl1_clamped = {
        spl1_clamped.insert_knots(&as_)
    };

    let spl1_clamped = {
        rebase_at_left(&spl1_clamped, 0.0, ks.as_slice())
    };
    let spl2_clamped = {
        rebase_at_left(&spl2_clamped, 0.0, ks.as_slice())
    };

    let cpts1 = spl1_clamped.control_points();
    let cpts2 = spl2_clamped.control_points();

    let mut cpts: Vec<T> = Vec::with_capacity(cpts1.len() + cpts2.len());
    for &cpt in cpts1.iter().rev() {
        cpts.push(cpt);
    }

    for (cpt, cpt2) in cpts[cpts1.len() - join_cont - 1..].iter_mut().zip(cpts2.iter()) {
        *cpt = cpt.lerp(0.5, *cpt2)
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
    where T: PointT,
          <T as PointAsVector>::Vector: Norm<NormType = f64>
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

    let exs = rebase_at_right(&s, s.front(), ks.as_slice());
    let newks: Vec<f64> = exs.knots().iter().cloned().chain(once(1.0 + delta)).collect();
    let newcpts = exs.control_points().iter().cloned().chain(once(*target)).collect();

    Bspline::new(newcpts, newks)
}

#[test]
fn it_works() {}