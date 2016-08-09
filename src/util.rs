use itertools::Itertools;

// its a pity rust does not have upper_bound
pub fn upper_bound<T>(v: &Vec<T>, x: T) -> usize
    where T: PartialOrd<T>
{
    let mut cnt: usize = v.len();
    let mut ub : usize = 0;
    while cnt > 0 {
        let step = cnt / 2;
        let mid = ub + step;
        if x.ge(&v[mid]) {
            ub = mid + 1;
            cnt -= step + 1;
        } else {
            cnt = step;
        }
    }
    ub
}

// use the itertools merge
pub fn merge<U>(orig: &Vec<U>, taus: &Vec<U>) -> Vec<U>
    where U: PartialOrd<U> + Clone
{
    orig.iter().cloned().merge(taus.iter().cloned()).collect()
}
