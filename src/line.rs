use tol::Tol;
use vectorspace::VectorSpace;
use curve::{Curve, FiniteCurve};
use std::f64::INFINITY;
use std::f64;
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
        Line{start:*p1, dir: *dir}
    }

    pub fn new_joining(p1: &Point, p2: &Point) -> Line<Point>
    {
        let dir = *p2 - *p1;
        Self::new(p1,&dir)
    }

    pub fn closest_points( &self,
                           l2: &Line<Point>) -> (Point, Point)

    {
        let l1 = self;
        let d1 = l1.direction();
        let d2 = l2.direction();
        let r =  l1.start_pt() - l2.start_pt();
        let b = d1.dot(&d2);
        let c = d1.dot(&r);
        let f = d2.dot(&r);
        let d = 1.0 - b * b;

        if d.small() {
            (l1.start_pt(), l2.start_pt()+l2.direction()*f)
        }else {
            let p1 = l1.start_pt() +  d1 * (d1.dot(&(d2 * f - r)) / d);
            let p2 = l2.start_pt() +  d2 * (d2.dot(&(r - d1 * c)) / d);
            (p1, p2)
        }
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
            l1.start_pt() + d1 * (d1.dot(&(d2*f-r))/d)
        }
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
        let candidates = [
            self.eval(self.a),
            self.eval(self.b),
            self.l.closest_point(p) ];
        let mut mindist = (*p - candidates[0]).len();
        let mut c = candidates[0];
        for j in 1..3 {
            let pc = (*p - candidates[j]).len();
            if pc < mindist {
                mindist = pc;
                c = candidates[j];
            }
        }
        return c;
    }

}
