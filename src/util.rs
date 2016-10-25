use itertools::Itertools;

// its a pity rust does not have upper_bound
pub fn upper_bound<T>(v: &[T], x: T) -> usize
    where T: PartialOrd<T>
{
    let mut cnt: usize = v.len();
    let mut ub: usize = 0;
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

macro_rules! uniq_ts {
    ($ts:expr) => {
          $ts.map(|&x| -> Param { Param(x) })
            .dedup()
            .map(|Param(y)| {y})
  }
}

pub fn merge<U>(orig: &[U], taus: &[U]) -> Vec<U>
    where U: PartialOrd<U> + Clone
{
    orig.iter().cloned().merge(taus.iter().cloned()).collect()
}
