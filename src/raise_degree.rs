use vectorspace::PointT;
use bspline::Bspline;
use curve::BlossomCurve;
use tol::Param;
use smat::clamp_ends;
use splinedata::SplineData;
use itertools::Itertools;
pub fn raise_degree<T: PointT>(spl: &Bspline<T>) -> Bspline<T> {
    let mut s = spl.clone();
    clamp_ends(&mut s);
    let ts = s.knots();
    let d = spl.degree() as usize;
    let tlen = ts.len();
    // Vec::with_capacity(2 * s.knots().len());
    let tsiter = ts.iter();
    let uniqts: Vec<f64> = uniq_ts!(ts[d..tlen - d].iter()).collect();
    let new_knots: Vec<f64> = tsiter.clone()
                                    .merge(&uniqts)
                                    .cloned()
                                    .collect();

    let num_new_knots = new_knots.len();
    let p: usize = (spl.degree() + 1) as usize;
    let num_new_cpts = num_new_knots - p - 1;
    let mut new_cpts: Vec<T> = Vec::with_capacity(num_new_cpts);
    for i in 0..num_new_cpts {
        let mut cv = T::splat(0.0);
        let ts = &new_knots[i + 1..i + 1 + p];
        for j in 0..p {
            let mut nts: Vec<f64> = Vec::with_capacity(p);
            nts.push(new_knots[i]);
            for l in 0..j {
                nts.push(ts[l]);
            }
            for l in (j + 1)..p {
                nts.push(ts[l]);
            }
            cv += spl.blossom_eval(0, &nts).to_vector();
        }
        cv *= 1.0 / p as f64;
        new_cpts.push(cv);
    }
    Bspline::new(new_cpts, new_knots)
}


#[test]
fn it_works() {
    use point::Pt1;
    use nalgebra::Norm;
    use tol::Tol;
    use curve::Curve;
    let spl = Bspline::new(vec![Pt1::new(0.0), Pt1::new(1.0), Pt1::new(0.0)],
                           vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
    let s = raise_degree(&spl);
    assert!(s.degree() == spl.degree() + 1);
    assert!((s.eval(0.5) - spl.eval(0.5)).norm().small());
}