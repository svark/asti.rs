use vectorspace::PointT;
use bspline::Bspline;
use splinedata::{SplineData, greville};
use curve::BlossomCurve;
use std::cmp::{max, min};
fn tube_around_spline<P>(spl: Bspline<P>) -> (Vec<Vec<f64>>, Vec<Vec<f64>>)
    where P: PointT
{
    let num_cpts1 = spl.control_points().len();

    let dim = P::dim();
    let mins: Vec<Vec<f64>> = vec![ Vec::with_capacity(num_cpts1);dim];
    let maxs: Vec<Vec<f64>> = vec![ Vec::with_capacity(num_cpts1);dim];

    for i in 0..num_cpts1 {
        let u = greville(spl1, i);
        let up = spl.eval(u);
        let lp = spl.control_points(i);
        for j in 0..dim {
            let minp = min(up[j], lp[j]);
            let maxp = max(up[j], lp[j]);
            mins[j].push(minp);
            maxs[j].push(maxp);
        }
    }
    (mins, maxs)
}

fn curve_intersect<P>(spl1: Bspline<P>, spl2: BSpline<P>)
    where P: PointT
{
    let num_cpts1 = spl1.control_points().len();
    let d1 = spl1.degree();
    let num_cpts2 = spl2.control_points().len();
    let (mins, maxs) = tube_around_spline(spl1);

    for j in 0..num_cpts2 {
        let v = greville(spl2, j);
    }
}