use vectorspace::VectorSpace;
use curve::{Curve, FiniteCurve};


#[derive(Debug)]
pub struct  LegendreForm<P:VectorSpace> {
    a : Vec<P>,
    s: f64,
    e: f64,
}

impl<P:VectorSpace> LegendreForm<P>
{
    pub fn new(a: Vec<P>, s: f64, e: f64) -> LegendreForm<P>
    {
        LegendreForm{a:a, s:s, e:e}
    }
}

impl <P:VectorSpace> Curve for LegendreForm<P>
{
    type T = P;

    fn eval(&self,u : f64) -> P {
       self.eval_derivative(u, 0)
    }
    
    fn eval_derivative(&self, u: f64, der_order: u32) -> P
    {
        let u = (u - self.s)/(self.e - self.s);
        let u = 2.0*u  - 1.0; // map (s,e) -> (-1,1)
        let mut p1 = u;
        let mut p0 = 1.0;
        let mut v = P::splat(0.0);
        let mut coeffs = vec![0.0;self.a.len()];
        coeffs[0] = p0;
        coeffs[1] = p1;
        for i in 2..self.a.len() {
            let id = i as f64;
            // farouki 2000, bernstein to legendre conversion.
            let pnu = p1 * ((2*i-1) as f64) * u - p0 * (id - 1.0);
            coeffs[i] = pnu / id;
            p0 = p1;
            p1 = pnu; 
        }
        
        // use the recursion : D[P_{n+1}(x)] = (2n + 1) P_n(x) +  D[P_{n-1}(x)]
        // or in general D^j[P_{n+1}(x)] = (2n + 1) D^{j-1}P_n(x) +  D^j[P_{n-1}(x)]  
        for j in 1..(der_order as usize)
        {
            let mut coeffs_der : Vec<f64> = coeffs.clone();
            coeffs_der[j-1] = 0.0;
            for i in j..self.a.len() {
                let id = i as f64;
                coeffs_der[i] = (2.0*id - 1.0) * coeffs[i] + coeffs_der[i-1];
            }
            coeffs = coeffs_der.iter().cloned().map(|x| {x*2.0/(self.e-self.s)}).collect(); // account for remapped domain using chain rule
        }
        for (&c,&a) in coeffs.iter().zip(self.a.iter())
        {
            v = v + a*c;
        }
        v
    }
}

impl<P:VectorSpace> FiniteCurve for LegendreForm<P>
{
  fn param_range(&self) -> (f64, f64)
  {
      (self.s,self.e)
  }
}

#[test]
fn it_works() {
    use tol::Tol;
    let lf = LegendreForm::new(vec!(0.0,1.0,2.0), 0., 1.);
    println!("{:?}", lf.eval(0.0));
    assert!((lf.eval(0.0) - 1.0).len().small());
    println!("{:?}", lf.eval(1.0));
    assert!((lf.eval(1.0) - 3.0).len().small());

}
