use bspline::Bspline;
use curve::Curve;
use splinedata::{SplineData, KnotManip, greville};

pub fn find_next_root(spl: &mut Bspline<f64> , prev: f64, tol : f64) -> Option<f64>
{
    let p = spl.degree() as usize;
    let mut k = {
        let mut k = 1;
        let t = spl.knots();
        let c = spl.control_points();
        let n = c.len();
        while k < n && t[k] < prev
        {
            k = k + 1;
        }
        k
    };
    while let Some(root) = {
        let t = spl.knots();
        let c = spl.control_points();
        let n = c.len();

        while k < n && c[k-1] * c[k] > 0.0 {
            k += 1;
        }
        if k < n {
            let root = greville(spl,k)
                - c[k]*(t[k+p] - t[k])/(p as f64* (c[k] - c[k-1]));
            Some(root)
        }else {
            None
        }
    } {
        if spl.eval(root) < tol {
            return Some(root)
        }
        if spl.mult(root) >= p {
            k = k + 1;
        }else {
            *spl = spl.insert_knot(root);
        }
    }
    None
}

#[test]
fn it_works()
{
    use class_invariant::ClassInvariant;
    use tol::PARAMRES;
    use tol::Param;
    let mut spl  = Bspline::new(vec![-0.5,0.1, 0.2,-0.1], vec![0.,0.,0.,0.5,1.,1.,1.]);
    assert!(spl.is_valid().is_ok() );
    let mut lb = 0.0;
    let vs = vec![2.5/12.0, 11.0/12.0];
    {
        let mut j  = 0;
        while let Some(nroot) = find_next_root(&mut spl, lb, PARAMRES) {
            println!("r:{}",nroot);
            lb = nroot + 1.5*PARAMRES;
            assert_eq!(Param(vs[j]), Param(nroot));
            j+=1;
        }
        println!("no more roots");

    }
}