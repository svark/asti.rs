// #![feature(cfg_target_feature)]
// #[cfg(target_feature = "sse2")]
// extern crate simd;
//
extern crate itertools;
extern crate la;
extern crate nalgebra;
extern crate either;

pub mod errcodes;
pub mod point;
pub mod tol;
#[macro_use]pub mod util;
pub mod combinations;
#[macro_use]pub mod vectorspace;
pub mod angle;
pub mod rmat;
pub mod splinedata;
pub mod curve;
pub mod class_invariant;
pub mod bspline;
pub mod periodic_spline;
pub mod rational_bspline;
pub mod bezier;
pub mod smat;
pub mod monomial_form;
pub mod line;
pub mod rootfinder;
pub mod conic;
pub mod legendre_form;
pub mod change_basis;
pub mod param_finder;
pub mod integrate;
pub mod tessellate;
pub mod spline_approx;
pub mod box_compute;
pub mod closest_pt_on_curve;
pub mod raise_degree;
pub mod rev;
pub mod reparametrize;
pub mod trim_extend_join;
pub mod split_curve;
pub mod spline_interp;
pub mod fair_curve;
#[test]
fn it_works() {

    use vectorspace::{PointT, Ops};
    let p = point::Pt2::new(1.0, 1.0);
    let q = point::Pt2::new(2.0, 0.0);

    let s = point::Pt2::new(3.0, 5.0);
    assert!(p.extract(0) >= 1.0);
    assert!(q.extract(0) >= 2.0);
    assert!(s.extract(0) >= 3.0);
    let r = p.lerp(0.5, q) * 2.0;
    assert!(r.extract(0) >= 3.0);
    {
        let p = point::Pt4::new(1.0, 1.0, 1.0, 1.0);
        let q = point::Pt4::new(2.0, 0.0, 0.0, 0.0);

        let s = point::Pt4::new(3.0f64, 4.0f64, 5.0f64, 6.0f64);
        assert!(p.extract(0) >= 1.0);
        assert!(q.extract(0) >= 2.0);
        assert!(s.extract(0) >= 3.0);
        let r = p.lerp(0.5, q) * 2.0;
        assert!(r.extract(0) >= 3.0);
    }
}
