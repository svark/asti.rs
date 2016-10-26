use bspline::Bspline;
use splinedata::SplineData;
use vectorspace::PointT;
pub fn reverse_curve<T: PointT>(spl: &Bspline<T>) -> Bspline<T> {
    let rcpts = spl.control_points().iter().rev().cloned().collect::<Vec<_>>();
    let rts = spl.knots().iter().rev().map(|&x| -x).collect::<Vec<_>>();
    Bspline::new(rcpts, rts)
}

#[test]
fn it_works() {
    use point::Pt1;
    use tol::Tol;
    use curve::{Curve, FiniteCurve};
    use nalgebra::Norm;
    let spl = Bspline::new(vec![Pt1::new(0.0), Pt1::new(1.0), Pt1::new(2.1)],
                           vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
    let s = reverse_curve(&spl);
    assert!((s.start_param() + spl.end_param()).small());
    assert!((s.end_param() + spl.start_param()).small());
    assert!((s.eval(-1.) - spl.eval(1.)).norm().small());
    assert!((s.eval(-0.0) - spl.eval(0.0)).norm().small());
    assert!((s.eval(-0.5) - spl.eval(0.5)).norm().small());
}
