use splinedata::{KnotManip, SplineData, greville};
use smat::clamp_unclamped_ends;
use tol::Tol;
use bspline::Bspline;
use vectorspace::PointT;
use std::f64;
use curve::Curve;
enum MinMaxBound {
    MinBound,
    MaxBound,
}

fn find_knot_at_bound<P>(bs: &Bspline<P>, ci: i32, mm: &MinMaxBound) -> (f64, bool)
    where P: PointT
{
    let pts = bs.control_points();
    let t = bs.knots();
    let d = bs.degree() as usize;

    // asssumes that the curve is pre-clamped at ends
    assert_eq!(bs.end_mult(), d + 1);
    assert_eq!(bs.start_mult(), d + 1);

    macro_rules! ex {
        ($e:expr) => ($e.extract(ci as usize))
    }

    let (mut j, minmax) = match mm {
        &MinMaxBound::MinBound => {
            let mut min = (0usize, pts.first().unwrap());
            for (i, p) in pts.iter().enumerate() {
                if ex!(min.1) >= ex![p] {
                    min = (i, p);
                }
            }
            min
        }
        &MinMaxBound::MaxBound => {
            let mut max = (0usize, pts.first().unwrap());
            for (i, p) in pts.iter().enumerate() {
                if ex!(max.1) <= ex![p] {
                    max = (i, p);
                }
            }
            max
        }
    };

    j = (j..pts.len())
            .into_iter()
            .take_while(|&k| (ex!(pts[k]) - ex!(minmax)).small())
            .last()
            .unwrap();

    let i = (0..j + 1)
                .rev()
                .into_iter()
                .take_while(|&k| (ex!(pts[k]) - ex!(minmax)).small())
                .last()
                .unwrap();

    let v = greville(bs, i);

    if i == 0 || i == pts.len() - 1 {
        return (v, true);
    }


    if j - i + 1 >= d || j + 1 == pts.len() {
        // curve passes through this maxima
        return (v, true);
    }

    // insert a knot u for which we would have pts[i] ==
    // pts[j+1] after insertion
    let pi = ex!(pts[i]);
    let pj = ex!(pts[j + 1]);
    let pi1 = ex!(pts[i - 1]);

    let ai = (pi - pj) * (t[i + d] - t[i + 1]) /
             ((pi - pi1) * (t[j + d + 1] - t[j + 1]) + (pi - pj) * (t[i + d] - t[i]));
    let u = t[i + d] - ai * (t[i + d] - t[i]);
    (u, false)
}

fn find_bound_by_insertion<P: PointT>(bs: &Bspline<P>, i: i32, mm: MinMaxBound) -> (f64, f64) {
    let cbs = clamp_unclamped_ends(bs);
    let spl = cbs.as_ref().unwrap_or(bs);

    let (mut u, found_bound) = find_knot_at_bound(spl, i, &mm);
    if !found_bound {
        let mut newspl = spl.insert_knot(u);
        while {
            let (u_, f) = find_knot_at_bound(&newspl, i, &mm);
            u = u_;
            !f
        } {
            newspl = newspl.insert_knot(u);
        }
    }
    (u, bs.eval(u).extract(i as usize))
}

// find param at minima of the given spline projected at coord index `i'
pub fn find_min_bound_by_insertion<P: PointT>(bs: &Bspline<P>, i: i32) -> (f64, f64) {
    find_bound_by_insertion(bs, i, MinMaxBound::MinBound)
}

// find param at maxima of the given spline projected at coord index `i'
pub fn find_max_bound_by_insertion<P: PointT>(bs: &Bspline<P>, i: i32) -> (f64, f64) {
    find_bound_by_insertion(bs, i, MinMaxBound::MaxBound)
}

pub fn tight_box<P: PointT>(bs: &Bspline<P>) -> (P, P) {
    let dim = P::dim();
    let mut min = P::zero_pt();
    let mut max = P::zero_pt();
    for i in 0..dim {
        min[i] = find_bound_by_insertion(bs, i as i32, MinMaxBound::MinBound).1;
        max[i] = find_bound_by_insertion(bs, i as i32, MinMaxBound::MaxBound).1;
    }
    (min, max)
}

pub fn compute_box<P>(bs: &Bspline<P>) -> (P, P)
    where P: PointT
{
    let num_cpts = bs.control_points().len();

    let dim = P::dim();
    let mut minp: P = P::splat(f64::MAX);
    let mut maxp: P = P::splat(f64::MIN);

    for i in 0..num_cpts {
        let lp = bs.control_points()[i];
        for j in 0..dim {
            let lpj = lp[j];
            if minp[j] > lpj {
                minp[j] = lpj;
            }
            if maxp[j] < lpj {
                maxp[j] = lpj;
            }
        }
    }
    (minp, maxp)
}

#[test]
fn it_works() {
    use bspline::Bspline;
    use point::Pt1;
    let ks = vec![0., 0., 0., 0., 0.0, 0.3, 0.6, 0.9, 1.0, 1.0, 1.0, 1.0, 1.0];
    let pts = vec![Pt1::new(0.0),
                   Pt1::new(0.1),
                   Pt1::new(0.45),
                   Pt1::new(0.7),
                   Pt1::new(0.66),
                   Pt1::new(0.53),
                   Pt1::new(0.41),
                   Pt1::new(0.3)];
    let bs = Bspline::new(pts, ks);
    let maxb = find_max_bound_by_insertion(&bs, 0);
    let minb = find_min_bound_by_insertion(&bs, 0);
    println!("param at maxima{:?}", maxb.0);
    println!("param at minima{:?}", minb.0);
    assert!((maxb.0 - 0.5786249416066778).abs().small());
    assert!((minb.0 - 0.0).abs().small());
    let (minb, maxb) = tight_box(&bs);

    assert_eq!(minb[0], 0.);
    assert_eq!(maxb[0], 0.6515263869973378);
}
