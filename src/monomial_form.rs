use curve::Curve;
use vectorspace::VectorSpace;
pub struct  MonomialForm<P: VectorSpace>
{
    pts: Vec<P>,
    s: f64,
    e: f64
}

impl<P:VectorSpace> MonomialForm<P>
{
    pub fn new(pts: Vec<P>, s: f64, e: f64) -> Self
    {
        MonomialForm{ pts: pts, s:s , e:e }
    }
}

impl<P:VectorSpace> Curve for MonomialForm<P>
{
    type T = P;
    fn param_range(&self) -> (f64,f64) {(self.s,self.e) }
    fn eval(&self,u : f64) -> P {
        let mut v : P = Default::default();
        let mut ui = 1.0;
        let u = (u - self.s)/(self.e - self.s);
        for i in 0..self.pts.len()
        {
            v = v +  self.pts[i] * ui;
            ui = ui*u;
        }
        v
    }

    fn eval_derivative(&self, u : f64, der_order: u32) -> P {
        let mut v : P = Default::default();
        let mut ui = 1.0;
        let u = (u - self.s)/(self.e - self.s);
        let mut f = 1 as usize;
        for j in 1..(der_order+1) as usize {
            f *= j;
        }
        let mut p = 1;
        for i in der_order as usize..self.pts.len()  {
            v  = v + self.pts[i] * (f as f64 * ui );
            ui *= u;
            f  *= (i+1)/p;
            p  += 1;
        }
        v
    }
}

#[test]
fn it_works()
{
    let mf = MonomialForm::new(vec![1.0;5], 0.0 ,1.0);
    assert_eq!(mf.eval(0.0) , 1.0);
    assert_eq!(mf.eval_derivative(0.0,2) , 2.0);
    assert_eq!(mf.eval_derivative(1.0,2) , 20.0);
}
