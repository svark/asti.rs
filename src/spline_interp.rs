pub use nalgebra::{Norm, Dot, Dimension, Absolute, Indexable};
use vectorspace::{PointT, to_pt, PV};
use la::Matrix;
use std::ops::Mul;
use tol::Tol;
use either::*;
use errcodes::GeomErrorCode;
use bspline::Bspline;
use periodic_spline::PeriodicBspline;
use bspline::SplineWrapper;
use std::iter::once;

#[derive(PartialEq)]
pub enum ParametrizationOption {
    ChordLength,
    CentripetalLength,
    AffinelyInvariant,
    NeilsonFoley,
}

#[derive(PartialEq)]
pub enum EndConditions {
    Periodic,
    ParabolicBlending,
    VanishingDoubleDerivatives,
    NotAKnot,
}

fn qmat<P>(pts: &[P]) -> Vec<f64>
    where P: PointT
{
    let mut n: usize = 0;
    let mut cpts = P::zero_pt();
    for q in pts {
        cpts.axpy(&1.0,q);
        n = n + 1;
    }
    let cg = cpts * (1.0 / (n as f64));
    let dim = P::dimension(None);
    let mut sigmaxy = vec![0.0;dim*dim];

    for &q in pts.into_iter() {
        let v: P::V = q - cg;
        for i in 0..dim {
            for j in i..dim {
                sigmaxy[i * dim + j] += v[i] * v[j];
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

macro_rules! psums {
  {$y:expr}  => {
      once(0.0).chain($y).scan(0.0, |s,e| { *s += e; Some(*s) }).collect::<Vec<f64>>()
  }
}

pub fn find_parameters<P:PointT>(pts: &[P], opt: ParametrizationOption) -> Vec<f64>
{
    match opt {
        ParametrizationOption::ChordLength => {
            psums![pts.windows(2)
              .map(|p| (p[1] - p[0]).norm())]
        }
        
        ParametrizationOption::CentripetalLength => {
            psums![pts.windows(2)
              .map(|p| (p[1] - p[0]).norm().sqrt())
            ]
        }
        
        ParametrizationOption::AffinelyInvariant => {
            let n = P::dimension(None);
            let ref m = Matrix::new(n, n, qmat(pts));
            psums![ pts.windows(2)
              .map(|p| to_pt::<P>(p[1] - p[0]))
              .map(|y| {
                  let ref maty = as_mat![y, n];
                  let my = m.mul(maty);
                  maty.t().dot(&my).sqrt()
              })]
        }

        ParametrizationOption::NeilsonFoley => {
            let n = P::dimension(None);
            let ref m = Matrix::new(n, n, qmat(pts));
            let widths = pts.windows(2)
                           .map(|p| to_pt::<P>(p[1] - p[0]))
                           .map(|y| {
                               let ref maty = as_mat![y, n];
                               let my = m.mul(maty);
                               maty.t().dot(&my).sqrt()
                           })
                           .collect::<Vec<f64>>();
            let thetas = pts.windows(2)
                           .map(|p| (p[1] - p[0]))
                           .zip(pts.windows(2).skip(1).map(|p| (p[1] - p[0])))
                           .map(|(u, v)| (u.dot(&v) / (u.norm() * v.norm())).acos())
                           .collect::<Vec<f64>>();
            psums![
                     once(widths[0]).chain((1..n - 1)
                     .into_iter()
                     .map(|i| {
                         widths[i] *
                        (1.0 + 3.0 * thetas[i] * widths[i - 1] / 2.0 * (widths[i - 1] + widths[i]) +
                        3.0 * thetas[i + 1] * widths[i + 1] / 2.0 * (widths[i] + widths[i + 1]))
                    }))
                    .chain(once(widths[n-1]))
                ]
                 
        }
    }
}

macro_rules! e{
    ($t:expr, $j:expr) => {$t[$j+1] - $t[$j]};
    ($t:expr, $j:expr, $step:expr) => {$t[$j+$step] - $t[$j]}
}

fn setrow<P: PointT>(lhs: &mut Matrix<f64>, i: usize, pv: P) {
    let dim = P::dimension(None);
    for j in 0..dim {
        let c: f64 = pv[j];
        lhs.get_mut_data()[(i * dim + j)] = c;
    }
}

pub fn setup_mat<P:PointT>(m: &mut Matrix<f64>,
                    p: &[P],
                    t: &[f64],
                    vecs: Vec<(usize,<P as PV>::V)>,
                    d : &mut Matrix<f64>
                    )
{
    // set up rows 1 through n - 2 of the matrix
    // as in hoschek  pg 88
    // eqn 3.18. contd.
    let n = t.len();
    // let dim = P::dimension(None);
    for j in 1..(n - 1) {
        let aj = e!(t, j);
        let aj1 = e!(t, j - 1);

        m.set(j, j - 1, aj);
        m.set(j, j, 2.0 * (aj + aj1));
        m.set(j, j + 1, aj1);

        let a = 3.0 * aj / aj1;
        let c = 3.0 * (aj1 / aj);

        let u = P::zero_pt() + e!(p, j - 1) * a + e!(p, j) * c;
        setrow(d, j, u);
    }

    // use wherever available, explicit information about derivatives.
    for (j, v) in vecs.into_iter() {
        if !v.norm().small() {
            m.set(j, j - 1, 0.0);
            m.set(j, j, 1.0);
            m.set(j, j + 1, 0.0);
            let pv = P::zero_pt() + v;
            setrow(d, j, pv);
        }
    }
}



fn eval_tangents<P>(
    pts :&[P],
    tb: &[f64],
    explicit_tgts : Vec<(usize,<P as PV>::V)>,
    ec : EndConditions
) -> Option< Vec<<P as PV>::V> >
where P:PointT
{
    let dim = P::dimension(None);
    let o = P::zero_pt();
    let n = pts.len();
    let ref mut m = Matrix::zero(n, n);
    let ref mut rhs: Matrix<f64> = Matrix::zero(n, 1);
    let b = match ec {
        EndConditions::Periodic => {
            assert!((pts[0] - pts[n - 1]).norm().small());

            m.set(0, n - 2, e!(tb, 0));
            m.set(0, 1, e!(tb, n - 2));
            m.set(0, 0, 2.0 * (e!(tb, 0) + e!(tb, n - 2)));

            let a0 = m.get(0, n - 2) / m.get(0, 1);
            // p0' = pn'
            m.set(n - 1, n - 1, -1.0);
            m.set(n - 1, 0, 1.0);

            setup_mat(m, pts, tb, explicit_tgts, rhs);
            let p = o + (e!(pts, 0)) * (3.0 * a0) + (e!(pts, n - 2)) * (3.0 / a0);
            setrow(rhs, 0, p);
            setrow(rhs, n - 1, o);
            m.solve(rhs)
        }

        EndConditions::ParabolicBlending => {
            m.set(0, 0, 1.);
            m.set(n - 1, n - 1, 1.);

            setup_mat(m, pts, tb, explicit_tgts, rhs);
            // see 3.19 expanded hoschek
            let p = o + e!(pts, 0) * 1.5 + e!(pts, 1) * (-0.5);
            setrow(rhs, 0, p);
            let q = o + e!(pts, n - 2) * 1.5 + e!(pts, n - 3) * (-0.5);
            setrow(rhs, n - 1, q);
            m.solve(rhs)
        }

        EndConditions::VanishingDoubleDerivatives => {
            // pg 88 of hoschek..  set up equation m p' = n.p for C^2
            // continuity of cubic polynomials at params. eq 3.16 we need to
            // solve for the tangent p' ( n of them) we have assumed that at
            // param 0 and 1 the double derivatives vanish (no curvature
            // condition)
            let a0 = e!(tb, 0);
            m.set(0, 0, a0 * 2.0);
            m.set(0, 1, a0);

            let a1 = e!(tb, n - 2);
            m.set(n - 1, n - 2, a1);
            m.set(n - 1, n - 1, 2.0 * a1);  // there is a mistake here in hoschek

            setup_mat(m, pts, tb, explicit_tgts, rhs);
            let p = o + e!(pts, 0) * 3.0;
            setrow(rhs, 0, p);
            let q = o + e!(pts, n - 2) * (-3.0);
            setrow(rhs, n - 1, q);

            m.solve(rhs)
        }
        EndConditions::NotAKnot =>
        {
            let a0 = e!(tb, 0);
            let a1 = e!(tb, 1);
            m.set(0, 0,(1./a0)* ( 1./a0 + 1./a1) );
            m.set(0, 1, (1./a0 + 1./a1)*(1.0/a0 + 1.0/a1));

            let b0 = 2.0/(a0 * a0 * a0);
            let b1 = 3.0/(a1 * a0 * a0);
            let b2 = 1.0/(a1 * a1 * a1);

            let h0 = e!(tb, n - 3);
            let h1 = e!(tb, n - 2);

            /* this gets messy so I let mathematica do the dirty work:
            second derivative at x for the cubic:
            Y[i_, x_] := 6*P[i + 1]*(1/(t[i + 1] - t[i])^2 - (2*(x - t[i]))/(t[i + 1] - t[i])^3) + 6*P[i]*((2*(x - t[i]))/(t[i + 1] - t[i])^3 - 1/(t[i + 1] - t[i])^2) + 2*Q[i + 1]*((3*(x - t[i]))/(t[i + 1] - t[i])^2 - 1/(t[i + 1] - t[i])) + 2*Q[i]*((3*(x - t[i]))/(t[i + 1] - t[i])^2 - 2/(t[i + 1] - t[i]))
            third derivative
            X[i] := Y[i_, x_] := 6*P[i + 1]*(1/(t[i + 1] - t[i])^2 - (2*(x - t[i]))/(t[i + 1] - t[i])^3) + 6*P[i]*((2*(x - t[i]))/(t[i + 1] - t[i])^3 - 1/(t[i + 1] - t[i])^2) + 2*Q[i + 1]*((3*(x - t[i]))/(t[i + 1] - t[i])^2 - 1/(t[i + 1] - t[i])) + 2*Q[i]*((3*(x - t[i]))/(t[i + 1] - t[i])^2 - 2/(t[i + 1] - t[i]))

            Not a knot condition both second and third derivatives are same at the knot t[N - 2] (last but one)
            eqn := Eliminate[ {Y[N - 3, t[N - 2]]==Y[N - 2, t[N - 2], X[N - 3]==X[N - 2]}, {Q[N - 3]}]

            eqn[[1, 1]]
            >> P[ - 1 +N] (3/(t[ - 3 +N] - t[ - 2 +N]) + 2/(t[ - 2 +N] - t[ - 1 +N]))

            Table[ Coefficient[eqn[[1, 2]], P[ - i +N]]//Simplify, {i, 2, 3}]

            >>{((t[ - 3 + N] - t[ - 1 + N])^2 (2 t[ - 3 + N] - 3 t[ - 2 + N] + t[ - 1 + N]))/((t[ - 3 + N] - t[ - 2 + N])^3 (t[ - 2 + N] - t[ - 1 + N])), \
            (t[ - 2 + N] - t[ - 1 + N])^2/(t[ - 3 + N] - t[ - 2 + N])^3}

            Table[ Coefficient[eqn[[1, 2]], Q[ - i +N]]//Simplify, {i, 1, 2}]

            >> {( - t[ - 3 + N] + t[ - 1 + N])/( t[ - 3 + N] - t[ - 2 + N]), - ((t[ - 3 + N] - t[ - 1 + N])^2/(t[ - 3 + N] - t[ - 2 + N])^2)}
             */

            m.set(n - 1, n - 1, (h0 + h1)/h1); //note the sign
            m.set(n - 1, n - 2, (h0 + h1)*(h0 + h1)/(h0*h0));

            setup_mat(m, pts, tb, explicit_tgts, rhs);
            let p = o + e!(pts, 0) * (b0 + b1)
                     + e!(pts, 1) * b2;
            setrow(rhs, 0, p);
            
            let mut q = o;
            let a = 3./h0 + 2./h1;
            let b = (h0 + h1)*(h0 + h1) * ( - 2.*h0 + h1 )/(h0*h0*h0*h1);
            let c = -h1*h1/(h0*h0*h0);
            q.axpy(&a, &pts[n-1] );
            q.axpy(&b, &pts[n-2]);
            q.axpy(&c, &pts[n-3]);
            setrow(rhs, n - 1, q);
            m.solve(rhs)
        }
    };
    if let Some(ext) = b {
        let mut vs: Vec<<P as PV>::V> = Vec::with_capacity(n);
        for i in 0..n {
            let mut v: P = P::zero_pt();
            for j in 0..dim {
                v[j] = ext.get(i, j);
            }
            vs.push(v.to_vector());
        }
        Some(vs)
    } else {
        None
    }
}


pub fn pchip_preconditions<P>(
    pts : &[P], //random access iter type
    end_conditions : EndConditions,
    vecs : &Vec< (usize, <P as PV>::V) >
    ) -> Result<(), GeomErrorCode> where P: PointT 
{
    let num_pts = pts.len();
    if num_pts <= 1{
        return Err(GeomErrorCode::NotEnoughPointsForInterp)
    }
    if vecs.last().unwrap().0 > num_pts  {
         return Err(GeomErrorCode::MismatchedArraySizes);
    }
    if end_conditions  == EndConditions::Periodic {
        if num_pts > 1 && !(*pts.first().unwrap() - *pts.last().unwrap()).norm().small() 
        {
            return Err(GeomErrorCode::InvalidPeriodicData);
        }else if num_pts <=2
        {
            return Err(GeomErrorCode::NotEnoughPointsForPeriodicInterp);
        }
    }
    
    for ps in pts.chunks(2).into_iter()
    {
        if (ps[1] - ps[0]).norm().small() 
        {  
            return Err(GeomErrorCode::DuplicatePointData);
        }   
   }
   Ok(())
}

fn pchip_open<P>(pts :&[P],
                 params: &[f64],
                  tgts: &[<P as PV>::V]
                  )
-> Bspline<P> where P:PointT 
{
    let n = pts.len();
    let mut cpts = vec![ P::zero_pt(); 2*n];
    let mut knots = vec![0.0;  2*n + 4];

    cpts[0] = pts[0];
    let v = tgts[0];
    cpts[1] = pts[0] +  v * (e!(params, 0)/3.0);
    for i in 0..4 {
       knots[i] = params[0];
    }

    for  i in 1..n - 1
    {
        let v = tgts[i];
        cpts[2*i] = pts[i]  + v *  (e!(params, (i - 1))  * (-1.0/3.0));
        cpts[2*i + 1] = pts[i]  + v * (e!(params, i)/3.0) ;
        knots[2*i + 2] = params[i];
        knots[2*i + 3] = params[i];
    }

    let v = tgts[n-1];
    cpts[2*n - 2] = pts[n - 1] + v * (-e!(params, n - 2)/3.0);
    cpts[2*n - 1] = pts[n - 1];

    for i in 0..4 {
        knots[2*n + i]= params[n - 1];
    }
    Bspline::new( cpts, knots)
}

fn pchip_closed<P:PointT>(pts :&[P],
                   params: &[f64],
                   tgts: &[<P as PV>::V]
                   )
-> PeriodicBspline<P> 
{
    let n = pts.len();
    let mut cpts = vec![P::zero_pt();2*n];
    let mut knots = vec![0.0;2*n + 4];
//pg 180, 188 hoschek, sisl s1379
    for i in 0..2{
        knots[i] = e!(params, n-2);
    }
   
    for i in 0..n {
        knots[2*i + 2] = params[i];
        knots[2*i + 3] = params[i];
    }

    knots[2*n + 2]  = params[n - 1] + e!(params, 0);
    knots[2*n + 3]  = knots[2*n + 3];

    for i in 0..n {
        cpts[2*i]     = pts[i] +  tgts[i] * (-e!(knots, 2*i + 1, 2) / 3.0);
        cpts[2*i + 1] = pts[i] +  tgts[i] * (e!(knots, 2*i + 2, 2) / 3.0);
    }

    assert!((pts[n - 1] - pts[0]).norm().small());
   
    PeriodicBspline::from_spline(Bspline::new(cpts, knots))
}

pub fn pchip<P:PointT>(pts :&[P],
         explicit_tgts : Vec<(usize,<P as PV>::V)>,
         ec : EndConditions, 
         po : ParametrizationOption) 
         -> Option<Either<Bspline<P>, PeriodicBspline<P> > >
{
    let params = find_parameters(pts, po);
    assert_eq!(params.len() , pts.len());
    let is_periodic = ec == EndConditions::Periodic;
    if let Some(vs) = eval_tangents(
                       pts, params.as_slice(),
                       explicit_tgts, ec) {
        let spl =  if is_periodic
            {
                Right(pchip_closed(pts, params.as_slice(), vs.as_slice()))
            }else {
                Left(pchip_open(pts, params.as_slice(), vs.as_slice()))
            };
        Some(spl)
    }else {
        None
    }
}

#[test]
pub fn it_works()
{
    use point::Pt1;
    use curve::Curve;
   let pts = vec![Pt1::new(1.6), Pt1::new(2.0),  Pt1::new(2.5) , Pt1::new(3.0), Pt1::new(3.1), Pt1::new(3.25), Pt1::new(3.1) ];
   let po = ParametrizationOption::CentripetalLength;
   let ec = EndConditions::NotAKnot;
   
   if let Some(Left(bs)) = pchip(&pts, vec![], ec, po) {
      assert_eq!(bs.eval(0.0), Pt1::new(1.6));
      assert_eq!(bs.eval(0.2), Pt1::new(1.727056090664379));
      assert_eq!(bs.eval(1.2), Pt1::new(2.3917073863769334));
   }else {
       assert_eq!(1,0);
   }
    
//     auto bs = geom::piecewise_cubic_hermite_interp(ps.begin(), ps.end(), opts,
//                                                        std::vector<double>(6, 0));
//     REQUIRE(bs.eval(0) == Approx(0.0));
//     REQUIRE(bs.eval(0.2) == Approx(0.802717));
//     REQUIRE(bs.eval(0.3806) == Approx(2));
//     REQUIRE(bs.eval(0.4) == Approx(2.08463));
//     REQUIRE(bs.eval(0.5) == Approx(2.36436));
//                                   REQUIRE(bs.eval(0.6) == Approx(2.56361));
//     //  REQUIRE(bs.eval(0.8957) == Approx( 3.2499998098 ));
//     REQUIRE(bs.eval(0.9) == Approx(3.24958));
//     REQUIRE(bs.eval(1.1) == Approx(2.95042));
    //
}
