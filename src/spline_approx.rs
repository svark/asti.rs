use bspline::Bspline;
use tol::PARAMRES;
use la::{Matrix};
use rmat::basis;
use point::Pt1;
pub fn cubic_approx1d(f : &Fn(f64)->f64, t: Vec<f64> ) -> Bspline<Pt1>
{
    let p = 3;
    let n = t.len() - p - 1;
    assert!(n>4); // place atleast one internal knot.

    let mut pts : Vec<Pt1>  = vec![Pt1::new(0.0);n];
    // pg 173 lyche
    //(@file :file-name "./lee_quasi.pdf" :to "./lee_quasi.pdf" :display "intro to quasi interpolants")
    for j in 2..n-2 {
        let t5 = vec!{ t[j + 1], (t[j + 1] + t[j + 2])/2.0,
                        t[j + 2], (t[j + 2] + t[j + 3])/2.0,
                        t[j + 3] - PARAMRES/2.0 }; // ./media/cubic_approx0.png

        let mut data = vec![0.0;25];
        for i in  0..5 {
            let ref b = basis(&t,(t5[i], None),3);
            for k in 0..4 {// ./media/cubic_approx.png
                if i > 1 {
                    data[5*i + k + 1] = b[k];
                }else {
                    data[5*i + k] = b[k]
                }
            }
        }
        let mat = Matrix::new(5, 5, data);
        let rhs= vec![f(t5[0]), f(t5[1]), f(t5[2]), f(t5[3]), f(t5[4])];
        let rhs_as_mat = Matrix::new(5,1,rhs);
        if let Some(res) = mat.solve(&rhs_as_mat) {
            let res = res.get_data();
            if j==2 {
                pts[j-2] = Pt1::new(res[0]);
                pts[j-1] = Pt1::new(res[1]);
            }
            pts[j] = Pt1::new(res[2]);
            if j == n-3 {
                pts[j+1] = Pt1::new(res[3]);
                pts[j+2] = Pt1::new(res[4]);
            }
        }
    }
    Bspline::new(pts, t)
}

#[test]
fn it_works() {
    use curve::Curve;
       use rootfinder::find_next_root;
       let f = |x| { x*x*x - x + 1.0}; // roots : 1.32471795724475
       let ref mut spl = cubic_approx1d(&f,vec![-2.0,-2.0,-2.0,-2.0,0.,2.,2.,2.,2.]);
       println!("{:?},{:?},{:?}", spl.eval(0.0),spl.eval(-1.324718), spl.eval(2.0));
       let root = find_next_root(spl,-2.0, 1e-8).unwrap();
       println!("{:?}",root);
       assert!((root + 1.32471795724475).abs() < 1e-8);

    {
        use nalgebra::Norm;
        use tol::Tol;
        let f = |x| { x*x + x + 1.0}; 
       let ref mut spl = cubic_approx1d(&f,vec![-1.0,-1.0,-1.0,-1.0,0.,1.,1.,1.,1.]);
       println!("{:?},{:?},{:?}", spl.eval(-1.0),spl.eval(0.0), spl.eval(1.0));
       assert!((spl.eval(-1.0) - Pt1::new(1.0)).norm().small());
       assert!((spl.eval(0.0) - Pt1::new(1.0)).norm().small());
       assert!((spl.eval(1.0) - Pt1::new(3.0)).norm().small());
       assert!((spl.eval(0.75f64.sqrt()-0.5) - Pt1::new(1.5)).norm().small());
    }
}