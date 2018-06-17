
pub mod create;
pub mod read;


/* calculate bin given an alignment covering [beg,end) (zero-based, half-closed-half-open) */
/* BAI specific code: use min_shift = 3, n_levels = 5.
pub fn reg2bin(beg: u32, end: u32) -> u32
{
    let end = end - 1;
    if beg>>14 == end>>14 { return ((1<<15)-1)/7 + (beg>>14); }
    if beg>>17 == end>>17 { return ((1<<12)-1)/7 + (beg>>17); }
    if beg>>20 == end>>20 { return ((1<<9)-1)/7 + (beg>>20); }
    if beg>>23 == end>>23 { return ((1<<6)-1)/7 + (beg>>23); }
    if beg>>26 == end>>26 { return ((1<<3)-1)/7 + (beg>>26); }
    return 0;
}
*/
/* calculate the list of bins that may overlap with region [beg,end) (zero-based) */


/* calculate bin given an alignment covering [beg,end) (zero-based, half-close-half-open) */
fn reg2bin(beg: u64, mut end: u64, min_shift: u32, n_levels: u32) -> u32
{
    let mut l = n_levels;
    let mut s = min_shift;
    let mut t = ((1<<n_levels*3) - 1) / 7;

    end -= 1;

    while l > 0 {
        if beg >> s == end >> s {
            return t + (beg >> s) as u32;
        }

        l -= 1;
        s += 3;
        t -= 1<<l*3;
    }

    return 0;
}


/// Calculate the list of bins that may overlap with region [beg,end) (zero-based)
// store result in `bins`
fn reg2bins(beg: u64, mut end: u64, min_shift: u32, depth: u32, bins: &mut Vec<u32>) {
    bins.clear();
    let mut l = 0;
    let mut t = 0;
    let mut s = min_shift + depth * 3;

    end -= 1;
    while l <= depth {
        let b = t + (beg >> s);
        let e = t + (end >> s);
        bins.extend(b as u32 ..  (e+1) as u32);

        s -= 3;
        t += 1 << l*3;
        l += 1;
    }
}

/*
pub fn reg2bins(beg: u32, end: u32, bins: &mut Vec<u16>)
{
    let end = end - 1;
    bins.clear();
    bins.push(0);

    for k in 1 + (beg>>26) .. 1 + (end>>26) {
        bins.push(k as u16);
    }

    for k in  9 + (beg>>23) ..  9 + (end>>23) {
        bins.push(k as u16);
    }

    for k in 73 + (beg>>20) .. 73 + (end>>20) {
        bins.push(k as u16);
    }

    for k in 585 + (beg>>17) .. 585 + (end>>17) {
        bins.push(k as u16);
    }

    for k in 4681 + (beg>>14) .. 4681 + (end>>14) {
        bins.push(k as u16);
    }
}
*/


#[inline]
fn hts_bin_bot(bin: u32, n_levels: u32) -> u32
{
    let mut l = 0;
    let mut b = bin;

    // Compute the level of the bin
    while b > 0 {
        l += 1;
        b = hts_bin_parent(b);
    }

    return (bin - hts_bin_first(l)) << (n_levels - l) * 3;
}

#[inline]
fn hts_bin_first(l: u32) -> u32 {
    (((1<<(((l)<<1) + (l))) - 1) / 7)
}

#[inline]
fn hts_bin_parent(l: u32) -> u32 {
    (((l) - 1) >> 3)
}

#[inline]
fn hts_bin_level(bin: u32) -> u32
{
    let mut l = 0;
    let mut b = bin;

    // Compute the level of the bin
    while b > 0 {
        l += 1;
        b = hts_bin_parent(b);
    }

    return l;
}