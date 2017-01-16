pub use nalgebra::{PointAsVector, Norm, Dot, Dimension, Absolute, Indexable};
use vectorspace::{PointT, to_pt, Ops};
use la::Matrix;
use std::ops::Mul;

pub enum ParametrizationOption {
    ChordLength,
    CentripetalLength,
    AffinelyInvariant,
    NeilsonFoley,
}

pub enum EndConditions {
    Periodic,
    ParabolicBlending,
    VanishingDoubleDerivatives,
    NotAKnot,
}

pub fn qmat<P>(pb: &[P]) -> Vec<f64>
    where P: PointT + Ops
{
    let mut n: usize = 0;
    let mut cpb = P::zero_pt();
    for q in pb {
        cpb += q.to_vector();
        n = n + 1;
    }
    let cg = cpb * (1.0 / (n as f64));
    let dim = P::dimension(None);
    let mut sigmaxy = vec![0.0;dim*dim];

    for &q in pb.into_iter() {
        let p: P = to_pt(q - cg);
        for i in 0..dim {
            for j in i..dim {
                sigmaxy[i * dim + j] += p.extract(i) * p.extract(j);
            }
        }
    }
    for i in 0..dim {
        for j in i..dim {
            sigmaxy[i * dim + j] *= 1.0 / (n as f64);
            sigmaxy[j * dim + i] *= 1.0 / (n as f64);
        }
    }
    sigmaxy
}

macro_rules! as_mat {
    {$y:expr,$n:expr}  => {
        {
            let mut vecy = vec![0.0; $n];
            for i in 0..$n {
                vecy[i] = $y[i]
            }
            Matrix::new($n,1,vecy)
        }
    }
}

pub fn find_parameters<P>(pb: &[P], opt: ParametrizationOption) -> Vec<f64>
    where P: PointT + Ops,
          <P as PointAsVector>::Vector: Dot<f64> + Norm<NormType = f64> + Clone
{
    match opt {
        ParametrizationOption::ChordLength => {
            pb.chunks(2)
              .map(|p| (p[1] - p[0]).norm())
              .scan(0.0, |x, y| {
                  *x += y;
                  Some(*x)
              })
              .collect::<Vec<f64>>()
        }
        ParametrizationOption::CentripetalLength => {
            pb.chunks(2)
              .map(|p| (p[1] - p[0]).norm().sqrt())
              .scan(0.0, |x, y| {
                  *x += y;
                  Some(*x)
              })
              .collect::<Vec<f64>>()
        }

        ParametrizationOption::AffinelyInvariant => {
            let n = P::dimension(None);
            let ref m = Matrix::new(n, n, qmat(pb));
            pb.chunks(2)
              .map(|p| to_pt::<P>(p[1] - p[0]))
              .scan(0.0, |x, y| {
                  let ref maty = as_mat![y, n];
                  let my = m.mul(maty);
                  *x += maty.t().dot(&my).sqrt();
                  Some(*x)
              })
              .collect::<Vec<f64>>()
        }

        ParametrizationOption::NeilsonFoley => {
            let n = P::dimension(None);
            let ref m = Matrix::new(n, n, qmat(pb));
            let widths = pb.chunks(2)
                           .map(|p| to_pt::<P>(p[1] - p[0]))
                           .map(|y| {
                               let ref maty = as_mat![y, n];
                               let my = m.mul(maty);
                               maty.t().dot(&my).sqrt()
                           })
                           .collect::<Vec<f64>>();
            let thetas = pb.chunks(2)
                           .map(|p| (p[1] - p[0]))
                           .zip(pb.chunks(2).skip(1).map(|p| (p[1] - p[0])))
                           .map(|(u, v)| (u.dot(&v) / (u.norm() * v.norm())).acos())
                           .collect::<Vec<f64>>();
            (1..n - 1)
                .into_iter()
                .map(|i| {
                    widths[i] *
                    (1.0 + 3.0 * thetas[i] * widths[i - 1] / 2.0 * (widths[i - 1] + widths[i]) +
                     3.0 * thetas[i + 1] * widths[i + 1] / 2.0 * (widths[i] + widths[i + 1]))
                })
                .scan(0.0, |x, y| {
                    *x = *x + y;
                    Some(*x)
                })
                .collect::<Vec<f64>>()
        }
    }
}
