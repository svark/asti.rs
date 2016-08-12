use vectorspace::VectorSpace;
use rmat::{locate_nu, sdiv};
use splinedata::SplineData;
use tol::PARAMRES;

pub struct Smat<'a> {
    a: f64,
    b: f64,
    knots: &'a [f64],
    deg: u32,
}

impl<'a> Smat<'a> {
    pub fn new(a: f64, b: f64, knots: &[f64], deg: u32) -> Smat {
        let nu = locate_nu(a, deg as usize, knots);
        Smat {
            a: a,
            b: b,
            knots: &knots[nu..],
            deg: deg,
        }
    }

    pub fn reval<T>(&self, cachea: &[T]) -> Vec<T>
        where T: VectorSpace
    {
        let size: usize = self.deg as usize + 1;
        let t = self.knots;
        let mut cachec = Vec::with_capacity(cachea.len());
        for j in 0..size {
            let mut cacheb: Vec<T> = cachea[0..size].to_owned();
            for sz in (j + 1..size).rev() {
                let lambda = sdiv(t[j + 1 - sz] - self.a, self.b - self.a);
                for i in 0..sz {
                    cacheb[i] = cacheb[i].lerp(lambda, cacheb[i + 1]);
                }
            }
            for sz in (1..j + 1).rev() {
                let mu = sdiv(t[sz] - self.a, self.b - self.a);
                for i in 0..sz {
                    cacheb[i] = cacheb[i].lerp(mu, cacheb[i + 1]);
                }
            }
            cachec.push(cacheb[0]);
        }
        cachec
    }

    pub fn seval<T>(&self, cachea: &[T]) -> Vec<T>
        where T: VectorSpace
    {
        let size: usize = self.deg as usize + 1;
        let t = &self.knots;
        let mut cachec = Vec::with_capacity(cachea.len());
        for i in 0..size {
            let mut cacheb: Vec<T> = cachea[0..size].to_owned();
            for sz in ((i + 1)..size).rev() {
                for j in 1..sz + 1 {
                    let lambda = sdiv(self.a - t[j - sz], t[j] - t[j - sz]);
                    cacheb[j - 1] = cacheb[j - 1].lerp(lambda, cacheb[j]);
                }
            }
            for sz in (1..i + 1).rev() {
                for j in 1..sz + 1 {
                    let mu = sdiv(self.b - t[j - sz], t[j] - t[j - sz]);
                    cacheb[j - 1] = cacheb[j - 1].lerp(mu, cacheb[j]);
                }
            }
            cachec.push(cacheb[0])
        }
        cachec
    }
}


pub fn rebase_at_left<T>(spl: &T, a: f64, us: &[f64]) -> T
    where T: SplineData
{
    let t = spl.knots();
    let deg = spl.degree() as usize;
    let cpts = spl.control_points();

    let nu = locate_nu(a, deg, t);

    debug_assert!(us.iter().all(|&u| u <= a));
    let b = t[nu + 1];
    let cpts = {
        let sm_t = Smat::new(a, b, t, spl.degree());
        // get cpts wrt bernstein basis
        sm_t.seval(&cpts[(nu - deg)..])
    };

    let t = {
        // knots updated to rebase at `a'
        let mut ks = t[(nu - deg)..].to_owned();
        for j in 0..(deg + 2) {
            ks[j] = us[j];
        }
        ks
    };

    let cpts = {
        // get cpts wrt bspline basis
        let sm_ks = Smat::new(a, b, t.as_slice(), spl.degree());
        sm_ks.reval(cpts.as_slice())
    };

    T::new(cpts, t)
}

pub fn rebase_at_right<T>(spl: &T, b: f64, us: &[f64]) -> T
    where T: SplineData
{
    let t = spl.knots();
    let deg = spl.degree() as usize;
    let cpts = spl.control_points();
    let nu = locate_nu(b - PARAMRES / 2.0, deg, t);
    let a = t[nu];

    debug_assert!(us.iter().all(|&u| u >= b));

    // get cpts wrt bernstein basis
    let cpts = {
        let sm_t = Smat::new(a, b, t, spl.degree());
        sm_t.seval(&cpts[..(nu + 1)])
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
        let sm_ks = Smat::new(a, b, t.as_slice(), spl.degree());
        sm_ks.reval(cpts.as_slice())
    };
    T::new(cpts, t)
}

pub fn clamp_at_right<T>(b: f64, spl: &T)
    where T: SplineData
{
    rebase_at_right(spl, b, &vec![b; spl.degree() as usize + 1]);
}

pub fn clamp_at_left<T>(b: f64, spl: &T)
    where T: SplineData
{
    rebase_at_left(spl, b, &vec![b; spl.degree() as usize + 1]);
}
