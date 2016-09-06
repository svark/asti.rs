use monomial_form::MonomialForm;
use splinedata::SplineData;
use curve::FiniteCurve;
use vectorspace::VectorSpace;
use legendre_form::LegendreForm;
use bezier::Bezier;
use class_invariant::ClassInvariant;


pub fn   to_monomial<P:VectorSpace>( bezf : &Bezier<P> ) -> MonomialForm<P>
{
    debug_assert!(bezf.is_valid().is_ok());
    let b = bezf.control_points();
    let sz = b.len();

    let mut monf = vec![P::splat(0.0); sz];
    let mut cni = 1;
   
    let n = sz - 1;
    let mut signi : i32 = -1;
    for i in 0..sz
    {
        signi  = -signi;

        if i != 0 {
            cni  = cni*(n - i + 1);
            cni  = cni / i;
        }
        let mut cil : usize = 1;

        let mut signl = -1;
        for l in 0.. (i+1)
        {
            signl = -signl;
            if l != 0 {
                cil = cil * (i - l + 1);
                cil = cil/l;
            }
            let mut co = (cni * cil) as f64;
            if signi * signl < 0 
            {
                co = -co;
            }
            monf[i] = b[l] * co;
        }
    }
    let (s,e) = bezf.param_range();
    MonomialForm::new(monf, s, e)
}

pub fn   to_bezier<P:VectorSpace>( mf : &MonomialForm<P>) -> Bezier<P>
{
    let mflen = mf.len();
    let mut cpts = vec![P::splat(0.0); mflen];
    let n = mflen - 1;
    for l in 0..mflen {
        let mut clk : usize = 1;
        let mut cnk : usize = 1;
        for k in 0..l+1 {
            if k != 0 {
                cnk *= n - k +1;
                cnk /= k;
                clk *= l - k + 1;
                clk /= k;
            }
            cpts[l] = cpts[l] +  mf.point(k) * ((clk as f64) / (cnk as f64));
        }
    }
    let (s,e) = mf.param_range();
    Bezier::new(cpts, s, e)
}

#[inline]
fn ncks(n: usize)-> Vec<usize> 
{
    let mut ncks_ : Vec<usize> = Vec::with_capacity(n+1);
    let mut val: usize  = 1;
    ncks_.push(1);
    for k in 1..n+1
    {
        val *= n - k + 1;
        val /= k;
        ncks_.push(val);
    }
    ncks_
}

// return C^{n+k}_k values for k in [0,j]
#[inline]
fn nkcks(n:usize, j:usize) ->Vec<usize>
{
    let mut nkcks_ : Vec<usize> = Vec::with_capacity(j+1);
    let  mut val = 1;
    nkcks_.push(1);
    for k in 1 ..j+1
    {
        val *= n + k;
        val /= k;
        nkcks_.push(val);
    }
    nkcks_
}


pub fn  to_legendre<P:VectorSpace>(bezf : &Bezier<P>) -> LegendreForm<P>
{
    let b   = bezf.control_points();
    let sz  = (bezf.degree() + 1) as usize;
    
    let n    = bezf.degree() as usize;

    let mut legf: Vec<P> = vec![P::splat(0.0); sz];

    let njcns : Vec<usize> = nkcks(n,n);
    for k in 0 .. n+1 {
        let nkjicji:Vec<usize> = nkcks(n-k, n);
        let kicis: Vec<usize>  = nkcks(k,n);
        let mut jsign: i32 =  -1;
        for j in 0 ..n+1
        {
            jsign *=- 1;
            let f = (2.0* (j as f64) + 1.0).sqrt()/((n + j + 1) as f64);
            let mut mjk:f64  = 0.0;
            let mut isign:i32 = -1;
            let mut jci:usize  = 1;
            for  i in 0..(j+1)  {
                isign *= -1;
                if i!= 0 {
                    jci  *= j - i + 1;
                    jci /= i;
                }
                let m = (jci * kicis[i] * nkjicji[j-i]) as f64;
                mjk += if isign > 0  {m} else {-m};
            }
            mjk /= njcns[j] as f64;
            mjk *= if jsign > 0 { f } else{ -f };
            let legfj =  legf[j] + b[k] * mjk;
            legf[j] = legfj;
        }
    }
  
    let (s,e) = bezf.param_range();
    LegendreForm::new(legf,s,e)
}

pub fn  from_legendre<P:VectorSpace>(lf: &LegendreForm<P>) -> Bezier<P>
{
    let n  =  lf.len() - 1;
    let mut ksign =  -1;
    let ncjs =  ncks(n);
    let mut cpts : Vec<P> = Vec::with_capacity(n+1);
    for k in 0..(n+1)
    {
        ksign = -ksign;
        let f   = (2.0*(k as f64) + 1.0).sqrt();

        let kcis  =   ncks(k) ;
        let nkcjs =   ncks(n-k) ;

        for  j in 0..(n+1)
        {
            let mut mjk   = 0.0;
            let mut isign = -1;
            use std::cmp::min;
            for  i in 0.. min(j+1 ,k+1)
            {
                isign = -isign;
                if n - k + i >= j 
                {
                    let kisq = kcis[i] * kcis[i];
                    let nkcji = nkcjs[j-i];
                    let m = (kisq * nkcji) as f64;
                    mjk += if isign > 0 { m } else { -m  };
                }
            }
            mjk *= if ksign >0 {f } else { -f};
            mjk /= ncjs[j] as f64;
            let cptsj = cpts[j] + lf.coeffs()[k]*mjk ;
            cpts[j] = cptsj;
        }
    }
    let (s,e) = lf.param_range();
    Bezier::new(cpts, s,e)
}

///#[test]
fn it_works() {
    use curve::Curve;
    use tol::Tol;

 //("monomial to bezier") 
    {
        let mut mon:Vec<f64> = Vec::with_capacity(5);
        for i in 0..5 {
            mon.push((i as f64) + 1.0);
        }
        let  mf = MonomialForm::new(mon.clone(), 0.0, 1.0);
        let bzf = to_bezier(&mf);
        assert!(bzf.eval(0.1).eqres(mf.eval(0.1)));
        let mf_dual  = to_monomial(&bzf);
        let cfs      = mf_dual.points();
        assert!(cfs.len()  == mon.len());
        for i in 0 ..cfs.len()
        {
            println!("c:{:?},m:{:?}",cfs[i], mon[i]);
            assert!(cfs[i].eqres(mon[i]));
        }
    }
    //("bezier  to monomial") 
    {
        let mut c:Vec<f64> = Vec::with_capacity(5);
        for i in 0..5 {
            c.push((i as f64) + 1.0) ;
        }

        let bf = Bezier::new(c.clone(), 1.0, 2.0);
        let bf_dual = to_bezier (&to_monomial(&bf));
        let  cfs = bf_dual.control_points();
        assert!(cfs.len() == c.len());
        for i in 0 ..cfs.len()
        {
            assert!(cfs[i] .eqres(c[i]));
        }
    }
    //("legendre to bezier") 
    {
        let mut legf:Vec<f64> = Vec::with_capacity(5);
        for i in 0..5 {
            legf.push((i as f64) + 1.0 );
        }
        let  lf = LegendreForm::new(legf.clone(), 0.0, 1.0);
        let bzf      = from_legendre(&lf);
        assert!(bzf.eval(0.1).eqres(lf.eval(0.1)));
        assert!(bzf.eval(0.9).eqres(lf.eval(0.9)));
        let lf_dual  = to_legendre(&bzf);
        let cfs      = lf_dual.coeffs();
        assert!(cfs.len()  == legf.len());
        for i in 0 ..cfs.len()
        {
            assert!(cfs[i].eqres(legf[i]));
        }
    }
    //("bezier to legendre") 
    {
        let mut c:Vec<f64> = Vec::with_capacity(5);
        for i in 0..5 {
            c.push((i as f64) + 1.0 );
        }
        let bf = Bezier::new(c.clone(), 1.0, 2.0);
        let  lf = to_legendre(&bf);
        let  bf_dual = from_legendre(&lf);
        assert!(bf.eval(1.1).eqres(lf.eval(1.1)));
        assert!(bf.eval(1.9).eqres(lf.eval(1.9)));
        let   cfs = bf_dual.control_points();
        assert!(cfs.len() == c.len());
        for i in 0 ..cfs.len()
        {
            assert!(cfs[i].eqres(c[i]));
        }
    }
}