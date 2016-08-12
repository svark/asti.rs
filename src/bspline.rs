use tol::Tol;
use vectorspace::VectorSpace;
use util::merge;
use splinedata::{SplineData, KnotManip};
use curve::{Curve, BlossomCurve};
use rmat::{eval, locate_nu};

pub struct Bspline<Point: VectorSpace> {
    control_points: Vec<Point>,
    knots: Vec<f64>,
    deg: u32,
}


impl<P> SplineData for Bspline<P>
    where P: VectorSpace
{
    type T = P;
    fn control_points(&self) -> &Vec<P> {
        &self.control_points
    }
    fn knots(&self) -> &Vec<f64> {
        &self.knots
    }
    fn degree(&self) -> u32 {
        self.deg
    }

    fn new(cpts: Vec<P>, ks: Vec<f64>) -> Bspline<P> {
        let d = ks.len() - cpts.len() - 1;
        Bspline {
            control_points: cpts,
            knots: ks,
            deg: d as u32,
        }
    }

}


impl<T> KnotManip for Bspline<T>
    where T: VectorSpace
{
    fn start_mult(&self) -> usize {
        let f = self.front();
        self.knots.iter().take_while(|&x| (*x - f).small_param()).count()
    }

    fn end_mult(&self) -> usize {
        let l = self.back();
        self.knots.iter().rev().take_while(|&x| (*x - l).small_param()).count()
    }

    fn front(&self) -> f64 {
        self.knots[self.deg as usize]
    }

    fn back(&self) -> f64 {
        self.knots[self.knots.len() - 1 - (self.deg as usize)]
    }


    fn locate_nu(&self, u: f64) -> usize {
        locate_nu(u, self.deg as usize, &self.knots)
    }

    fn mult(&self, u: f64) -> usize {
        let nu = self.locate_nu(u);
        self.knots[..nu + 1].iter().rev().take_while(|&x| (*x - u).small_param()).count()
    }

    fn rebase(&self, taus: Vec<f64>) -> Bspline<T> {
        let d = self.deg as usize;
        let ncpts = taus.len() - d - 1;
        let mut cpts: Vec<T> = Vec::with_capacity(ncpts);
        for i in 0..ncpts {
            cpts.push(self.blossom_eval(0, &taus[i..]));
        }
        Bspline::new(cpts, taus)
    }

    fn insert_knot(&self, tau: f64) -> Bspline<T> {
        let nu = self.locate_nu(tau);
        let mut newts: Vec<f64> = self.knots[0..nu + 1].to_owned();
        newts.push(tau);
        newts.extend(self.knots[nu + 1..].iter());
        let spl = self.rebase(newts);
        spl
    }

    fn insert_knots(&self, taus: &Vec<f64>) -> Bspline<T> {
        let merged = merge(&self.knots, taus);
        self.rebase(merged)
    }
}

pub trait ClassInvariant
{
    fn is_valid(&self) -> Result<bool, &str>;
}

impl<P> ClassInvariant for Bspline<P>
    where P: VectorSpace
{
    fn is_valid(&self) -> Result<bool, &str> {
        if self.knots.len() != self.control_points.len() + self.deg as usize + 1 {
            return Err("bad degree");
        }
        if self.deg < 1 {

            return Err("zero degree");
        }

        if self.control_points.len() < 1 {
            return Err("too few cpts");
        }

        let t = &self.knots;
        for j in 1..t.len() {
            let i = j - 1;
            if t[i] > t[j] {
                return Err("knots not sorted");
            }
        }
        let mut nextj: usize = 0;
        for j in 0..t.len() {
            if j < nextj {
                continue;
            }

            let m = self.mult(self.knots[j]);
            if m > self.deg as usize + 1 {
                return Err("too many dups in knots");
            }

            nextj = j + m;
        }
        return Ok(true);
    }
}

impl<P> BlossomCurve for Bspline<P>
    where P: VectorSpace
{
    type T = P;
    fn blossom_eval(&self, der_order: u32, us: &[f64]) -> P {
        let d = self.deg as usize;
        let nu = self.locate_nu(us[0]);
        let basis = eval(&self.knots, us, Some(nu), self.deg, der_order);
        let mut r: P = Default::default();
        let cpts = &self.control_points[nu - d..];
        for (&x, &y) in cpts.iter().zip(basis.iter()) {
            r = r + x * y;
        }
        r
    }
}


impl<P> Curve for Bspline<P>
    where P: VectorSpace
{
    type T = P;
    fn param_range(&self) -> (f64, f64) {
        let t = &self.knots;
        let d = self.deg as usize;
        let ncpts = t.len() - d - 1;
        (t[d], t[ncpts])
    }

    fn eval(&self, u: f64) -> P {
        let d = self.deg as usize;
        self.blossom_eval(0, &vec![u;d+1])
    }

    fn eval_derivative(&self, u: f64, order: u32) -> P {
        let d = self.deg as usize;
        self.blossom_eval(order, &vec![u;d+1])
    }
}


#[test]
fn it_works() {
    let bs = Bspline {
        control_points: vec![0.0, 1.0, 0.5],
        knots: vec![1.0, 1.0, 1.0, 2.0, 2.0, 2.0],
        deg: 2,
    };
    assert_eq!(bs.is_valid(), Ok(true));
    assert_eq!(bs.eval(1.0), 0.0);
    assert_eq!(bs.eval(1.5), 0.625);
    assert_eq!(bs.eval(2.0), 0.5);

    assert_eq!(bs.blossom_eval(0, &[1.0, 1.0, 1.0]), 0.0);
    assert_eq!(bs.blossom_eval(0, &[2.0, 2.0, 2.0]), 0.5);
    assert_eq!(bs.blossom_eval(0, &[1.0, 1.0, 2.0]), 1.0);

    assert_eq!(bs.eval_derivative(1.0, 1), 2.0);

    let spl = bs.insert_knot(1.1);
    assert_eq!(spl.is_valid(), Ok(true));
    assert_eq!(spl.eval(1.0), 0.0);
    assert_eq!(spl.eval(1.5), 0.625);
    assert_eq!(spl.eval(2.0), 0.5);
}
