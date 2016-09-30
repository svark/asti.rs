use tol::Tol;
use vectorspace::VectorSpace;
use curve::{Curve, FiniteCurve};
use std::f64::INFINITY;
use std::f64;
use angle::cross;
use point::Point3;
pub struct Line<P:VectorSpace>
{
    start: P,
    dir : P,
}

impl <P:VectorSpace> Line<P>
{
    fn start_pt(&self) -> P { self.start }
    fn direction(&self) -> P {self.dir }
}

pub struct LineSeg<P:VectorSpace>
{
    l : Line<P>,
    a : f64,
    b : f64,
}


impl<P:VectorSpace> Curve for LineSeg<P>
{
    type T = P;
    fn eval(&self, u: f64) -> P { self.l.start + self.l.dir * u }
    fn eval_derivative(&self, u: f64, der_order : u32) -> P {
        match der_order
        {
            0 => self.eval(u),
            1 => self.l.dir,
            _ => Default::default()
        }
    }
}

impl<P:VectorSpace> FiniteCurve for LineSeg<P>
{
    fn param_range(&self) -> (f64,f64) { (self.a,self.b) }
}

impl <Point:VectorSpace> Line<Point> {
    pub fn new(p1: &Point, dir: &Point) -> Line<Point>
    {
        debug_assert!(!dir.len().small());
        Line{start:*p1, dir: dir.normalize().unwrap()}
    }

    pub fn new_joining(p1: &Point, p2: &Point) -> Line<Point>
    {
        let dir = *p2 - *p1;
        Self::new(p1,&dir)
    }


    pub fn closest_point(&self, p:& Point) -> Point
    {
        let l = self;
        let p1 = l.start_pt() + l.direction() * (*p - l.start_pt()).dot(&l.direction());
        return p1;
    }

    pub fn intersect_with_line(&self, l2:&Self) -> Point
    {
        let l1 = self;
        let d1 = l1.direction();
        let d2 = l2.direction();
        let r = l1.start_pt() - l2.start_pt();

        let b = d1.dot(&d2);
        let f = d2.dot(&r);

        let d = 1.0 - b*b;

        if  d.abs().small()  {
            Point::splat(INFINITY)
        }else {
            let xpt = l1.start_pt() + d1 * (d1.dot(&(d2*f-r))/d);
            debug_assert!(((xpt - l2.start_pt()).dot(&d2).abs() - (xpt - l2.start_pt()).len()).small() );
            xpt
            
        }
    }
}

pub fn closest_points( l1: &Line<Point3>,
                       l2: &Line<Point3>) -> (Point3, Point3)

{
    let u = l1.direction();// these are unit dirs
    let v = l2.direction();

    let a = l1.start_pt();
    let c = l2.start_pt();
    let r = c - a;
    let p = cross(&u,&v);
    let plen = p.len();
    if plen.small() { // directions of the lines are parallel
        let cp : Point3= a - u*((a - c).dot(&u));
        (a, cp)
    }else {
        // reduce to planar case.
        let r = r - p *( r.dot(&p)/p.lensq() );
        let rlen = r.len();
        let cosgamma = r.dot(&u)/rlen;
        let cosalpha = u.dot(&v);
        let cosbeta = r.dot(&v)/rlen;
        //          +
        //         /a\
        //        /   \
        //       u     v 
        //      /      \
        //     /        \
        //    /          \
        //   +g____r____b+
        //  a             c
        let s =  rlen * (cosgamma + cosbeta * cosalpha) /(1.0 - cosalpha *  cosalpha);
        let t = rlen * cosbeta + s * cosalpha;
        let p1 = a +  u * s ;
        let p2 = c +  v * t;
        (p1, p2)
    }
}

impl <P:VectorSpace> LineSeg<P>
{
    pub fn new_joining(p1: &P, p2: &P) -> LineSeg<P>
    {
        let dir = *p2 - *p1;
        Self::new(Line::new(p1, &dir), 0.0, 1.0)
    }

    pub fn new(l : Line<P>, a: f64, b: f64) -> LineSeg<P>
    {
        LineSeg{ l: l, a:a, b:b}
    }

    pub fn closest_point(&self, p:& P) -> P
    {
        let (c,d) = (
            self.eval(self.a),
            self.eval(self.b)
        );
        let cp = self.l.closest_point(p);
        let proj_in_seg = (cp - c).dot(cp - d) > 0;
        if  proj_in_seg 
             return cp;

        if (*p - d).len() < (*p - c).len() {
             d
        }else {
             c   
        }
}
