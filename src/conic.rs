use tol::Tol;
use tol::RESABS;
use rational_bspline::{RationalBspline, rational_derivatives_from_derivatives};
use vectorspace::{PointT, Ops, to_pt};
use curve::{Curve, FiniteCurve};
use std::f64::consts::PI;
use errcodes::GeomErrorCode;
use errcodes::GeomErrorCode::*;
use line::Line;
use point::{is_collinear, is_coplanar, Pt2, Pt3, Vec3, Vec2};
use angle::{Angle, perp_in_plane};
use bspline::{Bspline, SplineWrapper};
use nalgebra::{Norm, Dot};

pub struct ConicArc<P: PointT> {
    p: [P::H; 3],
}

pub enum ConicType {
    Ellipse,
    Hyperbola,
    Parabola,
}

impl<P: PointT> ConicArc<P> {
    pub fn with_weights(p: &[P], w: &[f64]) -> Self {
        let pws = [p[0].hdim(1.0) * w[0], p[1].hdim(1.0) * w[1], p[2].hdim(1.0) * w[2]];
        ConicArc { p: pws }
    }

    pub fn new(p1: P::H, p2: P::H, p3: P::H) -> Self {
        ConicArc { p: [p1, p2, p3] }
    }

    pub fn weight(&self, i: usize) -> f64
        where P::H: Ops
    {
        self.p[i].extract(P::dim())
    }

    fn make_conic_arc_non_parallel(p: [P; 3], v: [P; 2]) -> Result<ConicArc<P>, GeomErrorCode>
     {
        let (l1, l2) = (Line::new(&p[0], &v[0]), Line::new(&p[2], &v[1]));
        let l1 = try!(l1.ok_or(GeomErrorCode::DegenerateOrSmallConic));
        let l2 = try!(l2.ok_or(GeomErrorCode::DegenerateOrSmallConic));

        let cpt = l1.intersect_with_line(&l2);
        let s = p[1]; // shoulder point;
        debug_assert!(!cpt[0].is_infinite());

        let (l3, l4) = (Line::new_joining(&p[0], &p[2]), Line::new_joining(&cpt, &s));
        let l3 = try!(l3.ok_or(GeomErrorCode::DegenerateOrSmallConic));
        let l4 = try!(l4.ok_or(GeomErrorCode::DegenerateOrSmallConic));
        let q = l3.intersect_with_line(&l4);

        let a = (q - p[0]).norm() / (q - p[2]).norm();
        let u = a / (a + 1.0);

        let p1 = cpt;
        let mut w = (1.0 - u) * (1.0 - u) * (s - p[0]).dot(&(p1 - s)) +
                    u * u * (s - p[2]).dot(&(p1 - s));
        w /= 2.0 * u * (1.0 - u) * (p1 - p[1]).norm_squared();

        let cpts = [p[0], p1, p[2]];
        let weights = [1.0, w, 1.0];
        Ok(ConicArc::with_weights(&cpts, &weights))
    }

    fn make_conic_arc_parallel(p: [P; 3], v: P) -> Result<ConicArc<P>, GeomErrorCode>
    {
        let s = p[1]; // shoulder point
        let (l1, l2) = (Line::new(&p[0], &p[2]), Line::new(&s, &v));
        let l1 = try!(l1.ok_or(GeomErrorCode::DegenerateOrSmallConic));
        let l2 = try!(l2.ok_or(GeomErrorCode::DegenerateOrSmallConic));
        let q = l1.intersect_with_line(&l2);

        let a = (q - p[0]).norm() / (q - p[2]).norm();
        let u = a / (1.0 + a);
        let b = ((1.0 - u) * (1.0 - u) + u * u) / (2.0 * u * (1.0 - u));
        let p1: P = to_pt((s - q) * b);
        let ref p1s = p1 - s;
        let w = ((1.0 - u) * (1.0 - u) * (s - p[0]).dot(p1s) + (s - p[2]).dot(p1s)) /
                (2.0 * u * (1.0 - u) * (p1 - p[1]).norm_squared());
        let cpts = [p[0], p1, p[2]];
        let weights = [1.0, w, 1.0];
        Ok(ConicArc::with_weights(&cpts, &weights))
    }

    fn conic_arc_preconditions(p: [Pt3; 3], v: [Vec3; 2]) -> Result<(), GeomErrorCode> {
        if v[0].norm().small() || v[1].norm().small() {
            Err(TangentVectorsTooSmall)
        } else if is_collinear(&p[0], &p[1], &p[2], RESABS) {
            Err(DegenerateOrSmallConic)
        } else if !is_coplanar(&p[0], &p[1], &p[2], &(p[0] + v[0]), RESABS) {
            Err(VectorsNotInPlaneOfPoints)
        } else if !is_coplanar(&p[0], &p[1], &p[2], &(p[2] + v[1]), RESABS) {
            Err(VectorsNotInPlaneOfPoints)
        } else {
            Ok(())
        }
    }

    pub fn typ(&self) -> ConicType {
        let w = self.weight(1).abs();
        if (w - 1.).small() {
            ConicType::Parabola
        } else if w < 1. {
            ConicType::Ellipse
        } else {
            ConicType::Hyperbola
        }
    }

    pub fn split_conic_at_shoulder(&self) -> (P::H, P::H, P::H) {
        let w = self.weight(1);
        assert!((self.weight(0) - 1.).abs().small());
        assert!((self.weight(2) - 1.).abs().small());
        let mut q1 = self.p[0].lerp(0.5, self.p[1]);
        let mut r1 = self.p[2].lerp(0.5, self.p[1]);
        // 7.42, pg 314 in the NURBS book
        let s = q1.lerp(0.5, r1).proj().unwrap();
        let w1 = ((1.0 + w) / 2.0).sqrt();
        q1 = q1 * (w1 / q1.extract(P::dim()));
        r1 = r1 * (w1 / r1.extract(P::dim()));
        (q1, s, r1)
    }
}

impl<P: PointT> Curve for ConicArc<P> {
    type T = <P::H as PointT>::L;
    fn eval(&self, u: f64) -> Self::T {
        self.eval_derivative(u, 0)
    }

    fn eval_derivative(&self, u: f64, der_order: u32) -> Self::T {
        let mut ders: Vec<P::H> = Vec::with_capacity(der_order as usize + 1);
        for i in 0..der_order + 1 {
            let der = match i {
                0 => {
                    let p01 = self.p[0].lerp(u, self.p[1]);
                    let p12 = self.p[1].lerp(u, self.p[2]);
                    p01.lerp(u, p12)
                }
                1 => {
                    let p01 = self.p[0].lerp(u, self.p[1]);
                    let p12 = self.p[1].lerp(u, self.p[2]);
                    let p01d: P::H = to_pt(self.p[1] - self.p[0]);
                    let p12d: P::H = to_pt(self.p[2] - self.p[1]);
                    p12 + (p01d.lerp(u, p12d) - p01)
                }
                2 => {
                    let dv: P::H = to_pt(self.p[1].lerp(0.5, self.p[2]) -
                                         self.p[0].lerp(0.5, self.p[1]));
                    dv * 2.0
                }
                _ => P::H::zero_pt(),
            };
            ders.push(der);
        }
        *rational_derivatives_from_derivatives(&ders).last().unwrap()
    }
}

impl<P: PointT> FiniteCurve for ConicArc<P> {
    fn param_range(&self) -> (f64, f64) {
        (0.0, 1.0)
    }
}

pub fn make_conic_arc(p: [Pt3; 3], v: [Vec3; 2]) -> Result<ConicArc<Pt3>, GeomErrorCode> {
    type ConicArc3 = ConicArc<Pt3>;
    try!(ConicArc3::conic_arc_preconditions(p, v));
    if (v[0].dot(&v[1]).abs() - v[0].norm() * v[1].norm()).small() {
        ConicArc::make_conic_arc_parallel(p, v[0].to_point())
    } else {
        ConicArc::make_conic_arc_non_parallel(p, [v[0].to_point(), v[1].to_point()])
    }
}

pub fn make_conic_arc_planar(p: [Pt2; 3], v: [Vec2; 2]) -> Result<ConicArc<Pt2>, GeomErrorCode> {
    if (v[0].dot(&v[1]).abs() - v[0].norm() * v[1].norm()).small() {
        ConicArc::make_conic_arc_parallel(p, v[0].to_point())
    } else {
        ConicArc::make_conic_arc_non_parallel(p, [v[0].to_point(), v[1].to_point()])
    }
}

pub fn make_circular_arc(p: [Pt3; 3]) -> Result<ConicArc<Pt3>, GeomErrorCode> {
    let p10 = p[0] - p[1];
    let p12 = p[2] - p[1];
    let p02 = p[2] - p[0];

    let angle = try!(p10.try_angle(&p12).ok_or(DegenerateCircle));
    let w = try!(p02.try_normalize(RESABS).ok_or(DegenerateCircle));
    let mut v = try!(perp_in_plane(w, p).ok_or(DegenerateCircle));
    let theta = PI - angle;
    if (p[1] - p[0]).dot(&v) < 0.0 {
        v = v * (-1.0);
    }
    let tgts = [v * theta.sin() + w * theta.sin(), v * (-theta.sin()) + w * theta.cos()];
    make_conic_arc(p, tgts)

}

pub fn make_rbspline_from_conic(arc: &ConicArc<Pt3>) -> Result<RationalBspline<Pt3>, GeomErrorCode> {
    let w = arc.weight(1);
    type ConicArc3 = ConicArc<Pt3>;
    type RationalBspline3 = RationalBspline<Pt3>;
    let dim = 3;
    if w < -1.0 {
        let ref mut rev = arc.p[1].clone();
        let nv = -arc.p[1].extract(dim);
        rev.replace(dim, nv);
        let cra = ConicArc3::new(arc.p[2].clone(), *rev, arc.p[0].clone());
        return make_rbspline_from_conic(&cra);
    }

    let projpts = [arc.p[0].ldim(), arc.p[1].ldim(), arc.p[2].ldim()];

    if let Some(alpha) = (projpts[1] - projpts[0]).try_angle(&(projpts[2] - projpts[1])) {
        if w >= 1.0 || (w > 0.0 && alpha > PI / 3.0) {
            // single segment parabola or hyperbola
            let ks = vec![0., 0., 0., 1., 1., 1.];
            let cpts = vec![arc.p[0].clone(), arc.p[1].clone(), arc.p[2].clone()];
            let spl = Bspline::new(cpts, ks);
            Ok(RationalBspline3::from_spline(spl))
        } else if w < 0.0 && alpha > PI / 2.0 {
            let ks = vec![0., 0., 0., 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1., 1., 1.];

            let (q0, s0, r0) = arc.split_conic_at_shoulder();

            let c1arc = ConicArc3::new(arc.p[0], q0, s0);
            let c2arc = ConicArc3::new(s0, r0, arc.p[2]);
            let (q1, s1, r1) = c1arc.split_conic_at_shoulder();
            let (q2, s2, r2) = c2arc.split_conic_at_shoulder();
            let cpts = vec![arc.p[0], q1, s1, r1, s0, q2, s2, r2, arc.p[2]];
            let spl = Bspline::new(cpts, ks);
            Ok(RationalBspline::from_spline(spl))
        } else {
            let (q, s, r) = arc.split_conic_at_shoulder();
            let ks = vec![0., 0., 0., 0.5, 0.5, 1., 1., 1.];
            let cpts = vec![arc.p[0], q, s, r, arc.p[2]];
            let spl = Bspline::new(cpts, ks);
            Ok(RationalBspline::from_spline(spl))
        }
    } else {
        Err(DegenerateOrSmallConic)
    }
}

#[test]
fn it_works() {
    // test ellipse
    {
        use point::Pt2;
        let pts = [Pt2::new(0.0, 0.0), Pt2::new(0.4, 0.3), Pt2::new(1.0, 0.0)];
        let vs = [Vec2::new(0., 1.), Vec2::new(0.4, 0.5)];

        if let Ok(arc) = make_conic_arc_planar(pts, vs) {
            let p0 = arc.eval(0.);
            assert!((p0 - pts[0]).norm().small());
            let p1 = arc.eval(1.0);
            assert!((p1 - pts[2]).norm().small());
            let p05 = arc.eval(0.5);
            assert!((p05 - Pt2::new(0.691566265, 0.478915663)).norm().small());
        } else {
            assert!(false);
        }
    }
}
