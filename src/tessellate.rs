use vectorspace::PointT;
use bspline::Bspline;
use splinedata::SplineData;
use bezier::split_into_bezier_patches;
use curve::{Curve, CurvePoint, FiniteCurve};
use nalgebra::{Norm, PointAsVector};
use std::ops::Add;

pub fn tessellate<P: PointT>(eps: f64, spl: &Bspline<P>) -> Vec<CurvePoint<P>>
    where <P as PointAsVector>::Vector: Norm<NormType = f64>
    + Add<<P as PointAsVector>::Vector, Output=<P as PointAsVector>::Vector>
{
    let d = spl.degree();
    let f = (d * (d - 1)) as f64;
    let mut patch_pts = Vec::with_capacity(16);
    for bpatch in split_into_bezier_patches(spl) {
        let mut maxnrm: f64 = 0.0;
        for bi in bpatch.control_points().as_slice().windows(3) {
            let v = (bi[0] - bi[1]) + (bi[2] - bi[1]);
            let nrm = v.norm();
            if maxnrm < nrm {
                maxnrm = nrm;
            }
        }
        let delta = (8.0 * eps * (f * maxnrm).recip()).sqrt();
        let s = bpatch.start_param();
        let e = bpatch.end_param();
        let n = delta.recip().ceil() as usize;
        let delta = (n as f64).recip();
        for i in 0..n {
            let p = s + (i as f64) * delta * (e - s);
            let cp = CurvePoint {
                pnt: bpatch.eval(p),
                t: p,
            };
            patch_pts.push(cp);
        }
    }
    let ep = spl.end_param();
    patch_pts.push(CurvePoint {
        pnt: spl.eval(ep),
        t: ep,
    });
    patch_pts
}

#[test]
fn it_works() {
    use point::Pt2;
    use bspline;
    use line::LineSeg;
    let p1 = Pt2::new(-1., 1.5);
    let p2 = Pt2::new(0., 1.);
    let p3 = Pt2::new(1., 0.);

    let bs = bspline::Bspline::new(vec![p1, p2, p3], vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
    let points: Vec<CurvePoint<Pt2>> = tessellate(1.0e-3, &bs);
    for p in points.as_slice().windows(2) {
        let t = 0.5 * (p[0].t + p[1].t);
        let ls: LineSeg<Pt2> = LineSeg::new_joining(&p[0].pnt, &p[1].pnt).unwrap();
        let cpt = bs.eval(t);
        let chord_pt = ls.closest_point(&cpt);
        assert!((cpt - chord_pt).norm() < 1.0e-3);
    }
}