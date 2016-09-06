use monomial_form::MonomialForm;
use change_basis::{to_bezier, to_monomial};
fn find_param_2d(p: &Point2, mf: MonomialForm<Point3>) -> Option<f64>
{
    let coeffs = Vec::with_capacity(mf.len());
    let a = |i: usize| { mf[i][0]};
    let b = |i: usize| { mf[i][1]};
    let c = |i: usize| { mf[i][2]};

    let l0j = Vec::with_capacity(mf.len());
    for m = 1..mf.size()
    {
        let j = m - 1;

        l0j.push( p[0] * (b(m) * c(0) - c(m) * b(0) ) + 
                 p[1] * ( a(0)* c(m) - c(0)*a(m)) + 
                 p[2] * ( a(m) * b(0) - b(m)*a(0)) );
    }
   let bzf = to_bezier(mf);
   Some(0.0)
}