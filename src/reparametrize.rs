use bspline::{Bspline, change_knots};
use vectorspace::PointT;
use splinedata::{SplineData, KnotManip};

pub fn reparametrize_start<T: PointT>(spl: Bspline<T>, a: f64) -> Bspline<T> {
    let prev = spl.front();
    let newknots: Vec<f64> = spl.knots().iter().cloned().map(|x| x - prev + a).collect();
    change_knots(newknots, spl)
}

pub fn reparametrize<T: PointT>(spl: Bspline<T>, a: f64, b: f64) -> Bspline<T> {
    let preva = spl.front();
    let prevb = spl.back();
    let newknots: Vec<f64> = spl.knots()
                                .iter()
                                .cloned()
                                .map(|x| a + (x - preva) * (b - a) / (prevb - preva))
                                .collect();
    change_knots(newknots, spl)
}

#[test]
fn it_works() {
    use point::Pt2;
    use tol::Tol;
    use nalgebra::Norm;
    use curve::Curve;

    let pts = vec![Pt2::new(0.0, 0.0),
                   Pt2::new(1.0, 1.0),
                   Pt2::new(1.5, 0.3),
                   Pt2::new(1.8, 0.1),
                   Pt2::new(2.0, 0.0)];
    let ks = vec![0.0, 0.0, 0.0, 0.5, 0.8, 1.0, 1.0, 1.0];
    let mut bs = Bspline::new(pts, ks.clone());
    bs = reparametrize_start(bs, 1.0);
    for (t, u) in ks.iter().cloned().zip(bs.knots().iter().cloned()) {
        assert!((u - t - 1.0).small());
    }
    assert!((bs.eval(1.0) - Pt2::new(0.0, 0.0)).norm().small());
    assert!((bs.eval(2.0) - Pt2::new(2.0, 0.0)).norm().small());

    bs = reparametrize(bs, 0.0, 1.0);
    for (t, u) in ks.iter().cloned().zip(bs.knots().iter().cloned()) {
        assert!((u - t).small());
    }
    assert!((bs.eval(0.0) - Pt2::new(0.0, 0.0)).norm().small());
    assert!((bs.eval(1.0) - Pt2::new(2.0, 0.0)).norm().small());
}
