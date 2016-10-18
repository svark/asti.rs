use curve::FiniteCurve;
use splinedata::{KnotManip, SplineData, greville};
use smat::clamp_ends;
use bspline::SplineWrapper;
use vectorspace::Ops;
use tol::{Tol, OSCoord};
enum MinMaxBound {
    MinBound,
    MaxBound,
}

fn find_knot_at_bound<SplineT>(bs: &SplineT, ci: i32, mm: &MinMaxBound) -> (f64, bool)
    where SplineT: SplineData + KnotManip,
          <SplineT as SplineData>::T: Ops
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

    let (j, minmax) = match mm {
        &MinMaxBound::MinBound => {
            pts.iter().enumerate().min_by_key(|&(_, p)| OSCoord(ex![p])).unwrap()
        }
        &MinMaxBound::MaxBound => {
            pts.iter().enumerate().max_by_key(|&(_, p)| OSCoord(ex![p])).unwrap()
        }
    };
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
    // pts[i+j] after insertion
    let pi = ex!(pts[i]);
    let pj = ex!(pts[j + 1]);
    let pi1 = ex!(pts[i - 1]);

    let ai = (pi - pj) * (t[i + d] - t[i + 1]) /
             ((pi - pi1) * (t[i + d + 1] - t[i + 1]) + (pi - pj) * (t[i + d] - t[i]));
    let u = t[i + d] - ai * (t[i + d] - t[i]);
    (u, false)
}

fn find_bound_by_insertion<T>(bs: &T, i: i32, mm: MinMaxBound) -> f64
    where T: SplineData + SplineWrapper + FiniteCurve + KnotManip + Clone,
          <T as SplineData>::T: Ops
{
    let mut spl = clamp_ends(bs.clone());
    let (mut u, mut found_bound) = find_knot_at_bound(&spl, i, &mm);
    while !found_bound {
        let newspl = spl.insert_knot(u);
        let (u_, f_) = find_knot_at_bound(&newspl, i, &mm);
        u = u_;
        found_bound = f_;
        spl = newspl;
    }
    bs.eval(u).extract(i as usize)
}

pub fn find_min_bound_by_insertion<T>(bs: &T, i: i32) -> f64
    where T: SplineData + SplineWrapper + FiniteCurve + KnotManip + Clone,
          <T as SplineData>::T: Ops
{
    find_bound_by_insertion(bs, i, MinMaxBound::MinBound)
}

pub fn find_max_bound_by_insertion<T>(bs: &T, i: i32) -> f64
    where T: SplineData + SplineWrapper + FiniteCurve + KnotManip + Clone,
          <T as SplineData>::T: Ops
{
    find_bound_by_insertion(bs, i, MinMaxBound::MaxBound)
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
    println!("param at{:?}", find_max_bound_by_insertion(&bs, 0));
    println!("param at minima{:?}", find_min_bound_by_insertion(&bs, 0));
    assert!((find_max_bound_by_insertion(&bs, 0) - 0.6515263864065358).abs().small());
    assert!((find_min_bound_by_insertion(&bs, 0) - 0.0).abs().small());
}