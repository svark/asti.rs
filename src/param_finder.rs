use bezier::Bezier;
use combinations::ncks;
use vectorspace::{PointT, Ops};
use point::{Pt2, Pt3};
use la::{Matrix, SVD};
use splinedata::SplineData;
use curve::{Domain};

// applications to cagd, thomas sederberg (book: cox,sturmfels,manocha, applns to computational algebraic geometry)
pub fn find_param_2d_weighted(p: &Pt2, bf: &Bezier<Pt3>) -> Option<f64> {

    let n = bf.degree() as usize;
    let mut l = vec![0.0;(n+1)*(n+1)];
    let (a, b) = (p.extract(0), p.extract(1));
    macro_rules! c {
        (x:$i:expr) => { bf.control_points()[$i].extract(0) };
        (y:$i:expr) => { bf.control_points()[$i].extract(1) };
        (w:$i:expr) => { bf.control_points()[$i].extract(2) };
    }

    let cn = ncks(n);
    for i in 0..(n + 1) {
        let cni = cn[i];
        let (xi, yi, wi) = (c!(x:i), c!(y:i), c!(w:i));
        for j in 0..i {
            let (xj, yj, wj) = (c!(x:j), c!(y:j), c!(w:j));
            let cnj = cn[j];
            let detx = a * (yi - yj) - b * (xi - xj) + (xi * yj - xj * yi);
            l[i * (n + 1) + j] = ((cni * cnj) as f64) * wi * wj * detx;
            l[j * (n + 1) + i] = -l[i * (n + 1) + j];
        }
    }

    let mut bigl = vec![0.0;n*n];
    for i in 0..n {
        for j in i..n {
            let mut biglij = 0.0;
            for k in 0..(i + 1) {
                let m = i + j + 1 - k;
                if m <= n {
                    biglij += l[k * (n + 1) + m];
                }
            }
            bigl[i * n + j] = biglij;
            bigl[j * n + i] = biglij;
        }
    }

    let m = Matrix::new(n, n, bigl);

    let svd = SVD::new(&m);
    let (sigma, v) = (svd.get_s(), svd.get_v());
    if sigma.get(n - 1, n - 1) > 1e-6 {
        return None;
    }
    let mut t = 0.0;
    let mut rs = 0.0;
    for j in 0..n {
        rs += v.get(j, n - 1);
    }
    let rsnf = (n as f64) * rs;
    for j in 0..n {
        t += v.get(j, n - 1) * (j as f64) / rsnf;
    }
    return Some(t);

}

pub fn find_param_2d(p: &Pt2, bf: &Bezier<Pt2>) -> Option<f64> {
    let lifted_pts = bf.control_points().iter().map(|&v| v.hdim(1.0)).collect::<Vec<_>>();
    find_param_2d_weighted(p,
                           &Bezier::new(lifted_pts, bf.param_range().0, bf.param_range().1))
}

#[test]
pub fn it_works() {
    use tol::Tol;
    use bspline::Bspline;
    use bezier::split_into_bezier_patches;
    use curve::Curve;
    use vectorspace::NVS;
    let spl = Bspline::new(vec![Pt3::new(0., 0., 1.0),
                                Pt3::new(0.0, 0.5, 1.0),
                                Pt3::new(0.5, 0.5, 1.0),
                                Pt3::new(0.5, 0., 1.0)],
                           vec![0., 0., 0., 0.5, 1., 1., 1.]);
    let pt_01 = spl.eval(0.1);
    println!("{:?}", pt_01);
    assert!(pt_01.extract(0) > 0.0);
    let p = Pt2::new(pt_01.extract(0), pt_01.extract(1));
    let mut found_t = false;
    for bp in split_into_bezier_patches(&spl).iter() {
        let q = bp.eval(0.1);
        assert!((pt_01 - q).norm().small());
        if let Some(t) = find_param_2d_weighted(&p, &bp) {
            assert!((t - 0.1).small());
            found_t = true;
        }
        break;
    }
    assert!(found_t);
}
