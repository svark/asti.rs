
#[inline]
pub fn ncks(n: usize) -> Vec<usize> {
    let mut ncks_: Vec<usize> = Vec::with_capacity(n + 1);
    let mut val: usize = 1;
    ncks_.push(1);
    for k in 1..n + 1 {
        val *= n - k + 1;
        val /= k;
        ncks_.push(val);
    }
    ncks_
}

// return C^{n+k}_k values for k in [0,j]
#[inline]
pub fn nkcks(n: usize, j: usize) -> Vec<usize> {
    let mut nkcks_: Vec<usize> = Vec::with_capacity(j + 1);
    let mut val = 1;
    nkcks_.push(1);
    for k in 1..j + 1 {
        val *= n + k;
        val /= k;
        nkcks_.push(val);
    }
    nkcks_
}
