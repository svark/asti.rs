use curve::{Curve, FiniteCurve};
use vectorspace::PointT;

#[derive(Debug)]
pub struct MonomialForm<P: PointT> {
    pts: Vec<P>,
    s: f64,
    e: f64,
}

impl<P: PointT> MonomialForm<P> {
    pub fn new(pts: Vec<P>, s: f64, e: f64) -> Self {
        MonomialForm {
            pts: pts,
            s: s,
            e: e,
        }
    }
    pub fn len(&self) -> usize {
        self.pts.len()
    }
    pub fn points(&self) -> &Vec<P> {
        &self.pts
    }

    pub fn point(&self, i: usize) -> P {
        self.pts[i]
    }
}

impl<P: PointT> Curve for MonomialForm<P> {
    type T = P;
    fn eval(&self, u: f64) -> P {
        let mut v: P = P::zero_pt();
        let mut ui = 1.0;
        let u = (u - self.s) / (self.e - self.s);
        for i in 0..self.pts.len() {
            v += (self.pts[i] * ui).to_vector();
            ui = ui * u;
        }
        v
    }

    fn eval_derivative(&self, u: f64, der_order: u32) -> P {
        let mut v: P = P::zero_pt();
        let mut ui = 1.0;
        let u = (u - self.s) / (self.e - self.s);
        let mut f = 1 as usize;
        for j in 1..(der_order + 1) as usize {
            f *= j;
        }
        let mut p = 1;
        for i in der_order as usize..self.pts.len() {
            v += (self.pts[i] * f as f64 * ui).to_vector();
            ui *= u;
            f *= i + 1;
            f /= p;
            p += 1;
        }
        v
    }
}

impl<P: PointT> FiniteCurve for MonomialForm<P> {
    fn param_range(&self) -> (f64, f64) {
        (self.s, self.e)
    }
}

#[test]
fn it_works() {
    use point::Pt1;
    let mf = MonomialForm::new(vec![Pt1::new(1.0);5], 0.0, 1.0);
    assert_eq!(mf.eval(0.0), Pt1::new(1.0));
    assert_eq!(mf.eval_derivative(0.0, 2), Pt1::new(2.0));
    assert_eq!(mf.eval_derivative(1.0, 2), Pt1::new(20.0));
}
