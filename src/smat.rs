use vectorspace::VectorSpace;
use rmat::{locate_nu, sdiv};

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
            deg: deg
        }
    }

    pub fn reval<T>(&self, cachea: &mut [T] )
        where T: VectorSpace
    {
        let size: usize = self.deg as usize + 1;
        let t = self.knots;
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
            cachea[j] = cacheb[0];
        }
    }

    pub fn seval<T>(&self, cachea: &mut [T])
        where T: VectorSpace
    {
        let size: usize = self.deg as usize + 1;
        let t = &self.knots;
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
            cachea[i] = cacheb[0]
        }
    }
}
