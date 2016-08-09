use point::lerp;
use bspline::{locate_nu, sdiv};
use std::ops::{Add, SubAssign, AddAssign, Mul};

pub struct Smat<'a> {
    a: f64,
    b: f64,
    knots: &'a [f64],
    deg: u32,
    offset: usize,
}

impl<'a> Smat<'a> {
    pub fn new(a: f64, b: f64, knots: &[f64], deg: u32) -> Smat {
        Smat {
            a: a,
            b: b,
            knots: knots,
            deg: deg,
            offset: locate_nu(a, deg as usize, knots),
        }
    }

    pub fn reval<T>(&self, cachea: &mut [T] )
        where T: Copy + Add<T, Output=T> + Mul<f64, Output = T> + AddAssign<T> + SubAssign<T> + Default
    {
        let nu = self.offset;
        let size: usize = self.deg as usize + 1;
        let t = &self.knots[nu..];
        let mut cachec: Vec<T> = Vec::new();
        cachec.resize(size, Default::default());
        for j in 0..size {
            let mut cacheb: Vec<T> = cachea[0..size].iter().cloned().collect();
            for sz in (j + 1..size).rev() {
                let lambda = sdiv(t[j + 1 - sz] - self.a, self.b - self.a);
                for i in 0..sz {
                    cacheb[i] = lerp(lambda, cacheb[i], cacheb[i + 1]);
                }
            }
            for sz in (1..j + 1).rev() {
                let mu = sdiv(t[sz] - self.a, self.b - self.a);
                for i in 0..sz {
                    cacheb[i] = lerp(mu, cacheb[i], cacheb[i + 1]);
                }
            }
            cachec[j] = cacheb[0];
        }
        for c in 0..cachec.len() {
            cachea[c] = cachec[c];
        }
    }

    pub fn seval<T>(&self, cachea: &mut [T])
        where T: Copy + Add<T, Output=T> + Mul<f64, Output = T> + AddAssign<T> + SubAssign<T> + Default
    {
        let nu = self.offset;
        let size: usize = self.deg as usize + 1;
        let t = &self.knots[nu..];
        let mut cachec: Vec<T> = Vec::new();
        cachec.resize(size, Default::default());
        for i in 0..size {
            let mut cacheb: Vec<T> = cachea[0..size].iter().cloned().collect();
            for sz in ((i + 1)..size).rev() {
                for j in 1..sz + 1 {
                    let lambda = sdiv(self.a - t[j - sz], t[j] - t[j - sz]);
                    cacheb[j - 1] = lerp(lambda, cacheb[j - 1], cacheb[j]);
                }
            }
            for sz in (1..i + 1).rev() {
                for j in 1..sz + 1 {
                    let mu = sdiv(self.b - t[j - sz], t[j] - t[j - sz]);
                    cacheb[j - 1] = lerp(mu, cacheb[j - 1], cacheb[j]);
                }
            }
            cachec[i] = cacheb[0]
        }
        for c in 0..cachec.len() {
            cachea[c] = cachec[c];
        }
    }
}
