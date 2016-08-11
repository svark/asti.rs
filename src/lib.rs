#![feature(cfg_target_feature)]
extern crate simd;
extern crate itertools;

pub mod point;
pub mod tol;
pub mod util;
pub mod vectorspace;
pub mod rmat;
pub mod splinedata;
pub mod curve;
pub mod bspline;
pub mod rational_bspline;
pub mod smat;
#[test]
fn it_works() {

    use vectorspace::VectorSpace;
    let p = point::Point2::new(1.0, 1.0);
    let q = point::Point2::new(2.0, 0.0);

    let b = vec![3.0f64, 4.0f64];
    let s = point::Point2::load(&b[..], 0);
    assert!(p.extract(0) >= 1.0);
    assert!(q.extract(0) >= 2.0);
    assert!(s.extract(0) >= 3.0);
    let r = p + q;
    assert!(r.extract(0) >= 3.0);
    {
        let p = point::Point4::new(1.0, 1.0, 1.0, 1.0);
        let q = point::Point4::new(2.0, 0.0, 0.0, 0.0);

        let b = vec![3.0f64, 4.0f64, 5.0f64, 6.0f64];
        let s = point::Point2::load(&b[..], 0);
        assert!(p.extract(0) >= 1.0);
        assert!(q.extract(0) >= 2.0);
        assert!(s.extract(0) >= 3.0);
        let r = p + q;
        assert!(r.extract(0) >= 3.0);
    }
}
