use vectorspace::PointT;
use curve::{Curve, FiniteCurve};


#[derive(Debug)]
pub struct LegendreForm<P: PointT> {
    a: Vec<P>,
    s: f64,
    e: f64,
}

impl<P: PointT> LegendreForm<P> {
    pub fn new(a: Vec<P>, s: f64, e: f64) -> LegendreForm<P> {
        LegendreForm { a: a, s: s, e: e }
    }

    pub fn coeffs(&self) -> &Vec<P> {
        &self.a
    }

    pub fn coeff(&self, k: usize) -> &P {
        &self.a[k]
    }

    pub fn len(&self) -> usize {
        self.a.len()
    }

    pub fn degree(&self) -> usize {
        self.a.len() - 1
    }
}

impl<P: PointT> Curve for LegendreForm<P> {
    type T = P;

    fn eval(&self, u: f64) -> P {
        self.eval_derivative(u, 0)
    }

    fn eval_derivative(&self, u: f64, der_order: u32) -> P {
        let u = (u - self.s) / (self.e - self.s);
        let u = 2.0 * u - 1.0; // map (s,e) -> (-1,1)
        let mut p1 = u;
        let mut p0 = 1.0;
        let mut v = P::splat(0.0);
        let mut coeffs = vec![0.0;self.a.len()];
        coeffs[0] = p0;
        coeffs[1] = p1;
        for i in 2..self.a.len() {
            let id = i as f64;
            // farouki 2000, bernstein to legendre conversion.
            let pnu = p1 * ((2 * i - 1) as f64) * u - p0 * (id - 1.0);
            coeffs[i] = pnu / id;
            p0 = p1;
            p1 = coeffs[i];
        }

        // use the recursion : D[P_{n+1}(x)] = (2n + 1) P_n(x) +  D[P_{n-1}(x)]
        // or in general D^j[P_{n+1}(x)] = (2n + 1) D^{j-1}P_n(x) +  D^j[P_{n-1}(x)]
        for j in 1..(der_order as usize) {
            let mut coeffs_der: Vec<f64> = coeffs.clone();
            coeffs_der[j - 1] = 0.0;
            for i in j..self.a.len() {
                let id = i as f64;
                coeffs_der[i] = (2.0 * id - 1.0) * coeffs[i] + coeffs_der[i - 1];
            }
            coeffs = coeffs_der.iter()
                               .cloned()
                               .map(|x| x * 2.0 / (self.e - self.s))
                               .collect(); // account for remapped domain using chain rule
        }
        for (i, (&a, &c)) in self.a.iter().zip(coeffs.iter()).enumerate() {
            v = v + (a * c * ((2 * i + 1) as f64).sqrt()).to_vector(); // sqrt(2i+1) makes legendre polynomials orthonormal
        }
        v
    }
}

impl<P: PointT> FiniteCurve for LegendreForm<P> {
    fn param_range(&self) -> (f64, f64) {
        (self.s, self.e)
    }
}

#[test]
fn it_works() {
    use tol::Tol;
    use point::Pt1;
    use nalgebra::Norm;
    let lf = LegendreForm::new(vec![Pt1::new(1.0), Pt1::new(0.0), Pt1::new(0.0)], 0., 1.);
    println!("{:?}", lf.eval(0.0));
    assert!((lf.eval(0.0) - Pt1::new(1.0)).norm().small());
    println!("{:?}", lf.eval(1.0));
    assert!((lf.eval(1.0) - Pt1::new(1.0)).norm().small());

    let lf = LegendreForm::new(vec![Pt1::new(0.0), Pt1::new(1.0), Pt1::new(0.0)], 0., 1.);
    println!("{:?}", lf.eval(0.0));
    assert!((lf.eval(0.0) - Pt1::new(-(3.0 as f64).sqrt())).norm().small());
    println!("{:?}", lf.eval(1.0));
    assert!((lf.eval(1.0) - Pt1::new((3.0 as f64).sqrt())).norm().small());

    let lf = LegendreForm::new(vec![Pt1::new(0.0), Pt1::new(0.0), Pt1::new(1.0)], 0., 1.);
    println!("{:?}", lf.eval(0.0));
    assert!((lf.eval(0.0) - Pt1::new((5.0 as f64).sqrt())).norm().small());
    println!("{:?}", lf.eval(1.0));
    assert!((lf.eval(1.0) - Pt1::new((5.0 as f64).sqrt())).norm().small());

}
