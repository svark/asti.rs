use vectorspace::PointT;
use rmat::{locate_nu, sdiv};
use splinedata::{SplineData, KnotManip};
use tol::PARAMRES;
use bspline::{Bspline, SplineWrapper};
use curve::FiniteCurve;

pub struct Smat<'a> {
    a: f64,
    b: f64,
    knots: &'a [f64],
    deg: u32,
    nu: usize,
}

impl<'a> Smat<'a> {
    pub fn new(a: f64, b: f64, knots: &[f64], nu: usize, deg: u32) -> Smat {
        Smat {
            a: a,
            b: b,
            knots: knots,
            deg: deg,
            nu: nu,
        }
    }

    pub fn reval<T>(&self, cachea: &[T]) -> Vec<T>
        where T: PointT
    {
        let size: usize = self.deg as usize + 1;
        let t = self.knots;
        let nu = self.nu;
        let mut cachec = Vec::with_capacity(cachea.len());
        for j in 0..size {
            let mut cacheb: Vec<T> = cachea[0..size].to_owned();
            for sz in (j + 1..size).rev() {
                let lambda = sdiv(t[j + 1 + nu - sz] - self.a, self.b - self.a);
                for i in 0..sz {
                    cacheb[i] = cacheb[i].lerp(lambda, cacheb[i + 1]);
                }
            }
            for sz in (1..j + 1).rev() {
                let mu = sdiv(t[sz + nu] - self.a, self.b - self.a);
                for i in 0..sz {
                    cacheb[i] = cacheb[i].lerp(mu, cacheb[i + 1]);
                }
            }
            cachec.push(cacheb[0]);
        }
        cachec
    }

    pub fn seval<T>(&self, cachea: &[T]) -> Vec<T>
        where T: PointT
    {
        let size: usize = self.deg as usize + 1;
        let t = &self.knots;
        let nu = self.nu;
        let mut cachec = Vec::with_capacity(cachea.len());
        for i in 0..size {
            let mut cacheb: Vec<T> = cachea[0..size].to_owned();
            for sz in ((i + 1)..size).rev() {
                for j in 1..sz + 1 {
                    let lambda = sdiv(self.a - t[j + nu - sz], t[j + nu] - t[j + nu - sz]);
                    cacheb[j - 1] = cacheb[j - 1].lerp(lambda, cacheb[j]);
                }
            }
            for sz in (1..i + 1).rev() {
                for j in 1..sz + 1 {
                    let mu = sdiv(self.b - t[j + nu - sz], t[j + nu] - t[j + nu - sz]);
                    cacheb[j - 1] = cacheb[j - 1].lerp(mu, cacheb[j]);
                }
            }
            cachec.push(cacheb[0])
        }
        cachec
    }
}


pub fn rebase_at_left<SplT>(spl: &SplT, a: f64, us: &[f64]) -> SplT
    where SplT: SplineWrapper
{
    let t = spl.knots();
    let deg = spl.degree() as usize;
    let cpts = spl.control_points();

    let nu = locate_nu(a + PARAMRES / 4.0, deg, t);

    debug_assert!(us.iter().all(|&u| u <= a));
    let b = t[nu + 1];
    let cpts = {
        let sm_t = Smat::new(a, b, t, nu, spl.degree());
        // get cpts wrt bernstein basis
        sm_t.seval(&cpts[(nu - deg)..])
    };

    let t = {
        // knots updated to rebase at `a'
        let mut ks = t[(nu - deg)..].to_owned();
        for j in 0..(deg + 1) {
            ks[j] = us[j];
        }
        ks
    };
    let mut cpts = {
        // get cpts wrt bspline basis
        let nu = locate_nu(a, deg, &t);
        let sm_ks = Smat::new(a, b, &t, nu, spl.degree());
        sm_ks.reval(&cpts)
    };
    if nu + 1 < spl.control_points().len() {
        cpts.extend(spl.control_points()[nu + 1..].to_owned());
    }
    SplT::from(Bspline::new(cpts, t))
}

pub fn rebase_at_right<T>(spl: &T, b: f64, us: &[f64]) -> T
    where T: SplineWrapper
{
    let t = spl.knots();
    let deg = spl.degree() as usize;
    let cpts = spl.control_points();
    let nu = locate_nu(b - PARAMRES / 4.0, deg, t);
    let a = t[nu];

    debug_assert!(us.iter().all(|&u| u >= b));

    // get cpts wrt bernstein basis
    let cpts = {
        let sm_t = Smat::new(a, b, t, nu, spl.degree());
        sm_t.seval(&cpts[nu - deg..])
    };
    // set up new knots
    let t = {
        let mut ks = t[..(nu + deg + 2)].to_owned();
        for j in 0..(deg + 1) {
            ks[nu + j + 1] = us[j];
        }
        ks
    };
    // get cpts wrt bspline basis
    let cpts = {
        let nu = locate_nu(a, deg, &t);
        let sm_ks = Smat::new(a, b, &t, nu, spl.degree());
        sm_ks.reval(cpts.as_slice())
    };
    if nu > deg {
        let mut lcpts = spl.control_points()[0..nu - deg].to_owned();
        lcpts.extend(cpts.into_iter());
        T::from(Bspline::new(lcpts, t))
    } else {
        T::from(Bspline::new(cpts, t))
    }
}

pub fn clamp_at_right<T>(b: f64, spl: &T) -> T
    where T: SplineWrapper
{
    rebase_at_right(spl, b, &vec![b; spl.degree() as usize + 1])
}

pub fn clamp_at_left<T>(a: f64, spl: &T) -> T
    where T: SplineWrapper
{
    rebase_at_left(spl, a, &vec![a; spl.degree() as usize + 1])
}

pub fn is_clamped<T>(spl: &T) -> bool
    where T: SplineData + KnotManip
{
    let d = spl.degree() as usize;
    spl.start_mult() == d + 1 && spl.end_mult() == d + 1
}

pub fn clamped_knots<T>(spl: &T) -> Vec<f64>
    where T: SplineData + KnotManip
{
    let d = spl.degree() as usize;
    let ts = spl.knots();
    let tlen = ts.len();
    let mut newts = vec![spl.front();d+1];
    newts.extend(&ts[d + 1..tlen - d - 1]);
    newts.extend(vec![spl.back();d+1]);
    newts
}

pub fn clamp_ends<T>(spl: T) -> T
    where T: SplineWrapper + FiniteCurve + KnotManip
{
    let mut clamped_spl = spl;
    let d = clamped_spl.degree() as usize;
    if clamped_spl.start_mult() != d + 1 {
        clamped_spl = clamp_at_left(clamped_spl.start_param(), &clamped_spl)
    }

    if clamped_spl.end_mult() != d + 1 {
        clamped_spl = clamp_at_right(clamped_spl.end_param(), &clamped_spl)
    }
    clamped_spl
}

pub fn clamp_unclamped_ends<T>(spl: &T) -> Option<T>
    where T: SplineWrapper + FiniteCurve + KnotManip
{
    let d = spl.degree() as usize;

    if spl.start_mult() != d + 1 {
        let clamped_spl = clamp_at_left(spl.start_param(), spl);
        Some(clamp_ends(clamped_spl))
    } else if spl.end_mult() != d + 1 {
        Some(clamp_at_right(spl.end_param(), spl))
    } else {
        None
    }
}


#[test]
fn it_works() {
    use bspline::Bspline;
    use splinedata::SplineData;
    use vectorspace::NVS;
    use curve::Curve;
    use point::Pt1;
    use tol::RESABS;
    let spl = Bspline::new(vec![Pt1::new(0.), Pt1::new(0.5), Pt1::new(0.)],
                           vec![0., 0., 0., 1., 1., 1.]);
    let pt_01 = spl.eval(0.1);
    let spl = rebase_at_left(&spl, 0., &vec![-0.2, -0.1, 0.0]);
    assert!((spl.eval(0.1) - pt_01).norm() < RESABS);
    let spl = rebase_at_left(&spl, 0., &vec![0., 0., 0.]);
    assert!(spl.control_points()
               .iter()
               .zip(vec![0., 0.5, 0.].iter())
               .all(|(&x, &y)| (x[0] - y).abs() < RESABS));
}
