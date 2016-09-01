use bezier::BezierForm;
use monomial_form::MonomialForm;
use splinedata::SplineData;
use curve::FiniteCuve;
use vectorspace::VectorSpace;

fn  <P:VectorSpace> to_monomial( bezf : BezierForm<P> ) -> MonomialForm<P>
{
    debug_assert!(bezf.is_valid().is_ok());
    let b = bezf.control_points();
    let sz = b.len();

    let mut monf = vec![P::splat(0.0); sz];
    let mut cni = 1;
   
    let n = sz - 1;
    for i 0..sz
    {
        let signi : i32 = (i & 1)<<1 
        let signi = signi - 1;

        if i != 0 {
            cni  = cni*(n - i + 1);
            cni  = cni / i;
        }
        let mut cil : usize = 1;

        for l in 0.. (i+1)
        {
            let signl : i32 = ((i & 1) << 1) - 1;
            if l != 0 {
                cil = cil * (i - l + 1);
                cil = cil/l;
            }
            let ev = &b[l];
            let mut co = cni * cil as f64;
            co  = co  * (signi * signl as f64);
            monf[i] = ev * co;
        }
    }
    let (s,e) = bezf.param_range();
    MonomialForm::new(monf, s, e)
}

fn <P:VectorSpace> to_bezier( mf : MonomialForm<P>) -> BeizerForm<P>
{
    let mflen = mf.len();
    let mut cpts = vec![P::splat(0.0); mflen]
    let n = mflen - 1;
    for l = 0..mflen {
        let clk : usize = 1;
        let cnk : usize = 1;
        for k = 0..l+1 {
            if k != 0 {
                cnk *= (n - k +1);
                cnk /= k;
                clk *- (l - k + 1);
                clk /= k;
            }
            cpts[l] = cpts[l] + (clk as f64) * mf.points(k) / (cnk as f64);
        }
    }
    let (s,e) = mf.param_range();
    BeizerForm::from_spline(Bspline::new(cpts, s, e);
}