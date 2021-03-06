use tol::Tol;
use vectorspace::{PointT, to_pt, PV, NVS, Cross};
use curve::{Curve, Domain};
use std::f64::INFINITY;
use std::f64;
use point::Pt3;
pub struct Line<P: PointT> {
    start: P,
    dir: P,
}

impl<P: PointT> Line<P> {
    fn start_pt(&self) -> &P {
        &self.start
    }
    fn direction(&self) -> <P as PV>::V {
        self.dir.to_vector()
    }
}

pub struct LineSeg<P: PointT> {
    l: Line<P>,
    a: f64,
    b: f64,
}


impl<P: PointT> Curve for LineSeg<P> {
    type T = P;
    fn eval(&self, u: f64) -> P {
        let &sp = self.l.start_pt();
        let ref d = self.l.dir;
        let mut res = sp.clone();
        res.axpy(&u, d);
        res
    }

    fn eval_derivative(&self, u: f64, der_order: u32) -> P {
        match der_order {
            0 => self.eval(u),
            1 => self.l.dir,
            _ => P::zero_pt(),
        }
    }
}

impl<P: PointT> Domain for LineSeg<P> {
    fn param_range(&self) -> (f64, f64) {
        (self.a, self.b)
    }
}

impl<Point: PointT> Line<Point> {
    pub fn new(p1: &Point, dir_as_pt: &Point) -> Option<Line<Point>> {
        let dir = (*dir_as_pt).as_vector();
        debug_assert!(!dir.norm().small());
        if let Some(nrml) = dir.try_normalize() {
            Some(Line {
                start: *p1,
                dir: Point::from_vec(&nrml),
            })
        } else {
            None
        }
    }

    pub fn new_joining(p1: &Point, p2: &Point) -> Option<Line<Point>> {
        let dir = to_pt(*p2 - *p1);
        Self::new(p1, &dir)
    }


    pub fn closest_point(&self, p: &Point) -> Point {
        let l = self;
        let d = l.direction();
        let f = (*p - l.start).dot(&d);
        let p1 = l.start + d * f;
        return p1;
    }

    pub fn intersect_with_line(&self, l2: &Self) -> Point {
        let l1 = self;
        let d1 = l1.direction();
        let d2 = l2.direction();
        let r = l1.start - l2.start;

        let b = d1.dot(&d2);
        let f = d2.dot(&r);

        let d = 1.0 - b * b;

        if d.abs().small() {
            Point::splat(INFINITY)
        } else {
            let f = d1.dot(&(d2 * f - r)) / d;
            let xpt = l1.start + d1 * f;
            xpt
        }
    }
}

pub fn closest_points(l1: &Line<Pt3>, l2: &Line<Pt3>) -> (Pt3, Pt3) {
    let u = l1.direction();// these are unit dirs
    let v = l2.direction();

    let &a = l1.start_pt();
    let &c = l2.start_pt();
    let r = c - a;
    let p = u.cross(&v);
    let plen = p.norm();
    if plen.small() {
        // directions of the lines are parallel
        let cp: Pt3 = a - u * ((a - c).dot(&u));
        (a, cp)
    } else {
        // reduce to planar case.
        let r = r - p * (r.dot(&p) / p.norm_squared());
        let rlen = r.norm();
        let cosgamma = r.dot(&u) / rlen;
        let cosalpha = u.dot(&v);
        let cosbeta = r.dot(&v) / rlen;
        //          +
        //         /a\
        //        /   \
        //       u     v
        //      /      \
        //     /        \
        //    /          \
        //   +g____r____b+
        //  a             c
        let s = rlen * (cosgamma + cosbeta * cosalpha) / (1.0 - cosalpha * cosalpha);
        let t = rlen * cosbeta + s * cosalpha;
        let p1 = a + u * s;
        let p2 = c + v * t;
        (p1, p2)
    }
}

impl<P: PointT> LineSeg<P> {
    pub fn new_joining(p1: &P, p2: &P) -> Option<LineSeg<P>> {
        let dir = to_pt(*p2 - *p1);
        if let Some(l) = Line::new(p1, &dir) {
            Some(Self::new(l, 0.0, 1.0))
        } else {
            None
        }
    }

    pub fn new(l: Line<P>, a: f64, b: f64) -> LineSeg<P> {
        LineSeg { l: l, a: a, b: b }
    }

    pub fn closest_point(&self, p: &P) -> P {
        let (c, d) = (self.eval(self.a), self.eval(self.b));
        let cp = self.l.closest_point(p);
        let proj_in_seg = (cp - c).dot(&(cp - d)) < 0.0;
        if proj_in_seg {
            return cp;
        }

        if (*p - d).norm() < (*p - c).norm() {
            d
        } else {
            c
        }
    }
}

#[test]
fn it_works() {
    use point::Pt2;
    let l = LineSeg::new_joining(&Pt2::new(1.0, 0.0), &Pt2::new(0.0, 0.0)).unwrap();
    assert_eq!(l.closest_point(&Pt2::new(0.5, 0.01)), Pt2::new(0.5, 0.0));
    assert_eq!(l.closest_point(&Pt2::new(1.1, 0.1)), Pt2::new(1.0, 0.0));
}