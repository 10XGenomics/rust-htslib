use nom::{le_u8, le_u32, le_i32, le_u64};
use std::path::Path;
use std::fs::File;
use std::io::Read;

use std::collections::HashSet;

use bam::IndexedReaderError;

#[derive(Copy, Clone, PartialEq, Eq, Ord, PartialOrd)]
pub struct Chunk {
    chunk_beg: u64,
    chunk_end: u64
}

struct Bin {
    bin: u32,
    chunks: Vec<Chunk>,
}

struct RefIndex {
    bins: Vec<Bin>,
    ioffsets: Vec<u64>,
}

struct BaiIndex {
    ref_indexes: Vec<RefIndex>,
    n_no_coor: u64,
}

named!(parse_ref_index<&[u8], RefIndex>,
    do_parse!(
        n_bin: le_i32 >>
        bins: count!(parse_bin, n_bin as usize) >>
        n_intv: le_i32 >>
        ioffsets: count!(le_u64, n_intv as usize) >>
        (
            RefIndex {
                bins: bins,
                ioffsets: ioffsets,
            }
        )
    )
);

named!(parse_chunk<&[u8], Chunk>,
  do_parse!(
    beg:          le_u64                 >>
    end:          le_u64                 >>
    (
      Chunk {
        chunk_beg: beg,
        chunk_end: end,
      }
    )
  )
);

named!(parse_bin<&[u8], Bin>,
    do_parse!(
        bin: le_u32 >>
        n_chunk: le_i32 >>
        chunks: count!(parse_chunk, n_chunk as usize) >>
        (
            Bin {
                bin: bin,
                chunks: chunks,
            }
        )
    )
);

const BAI_MAGIC: &[u8; 4] = &[b'B', b'A', b'I', 1u8];

named!(parse_bai<&[u8], BaiIndex>,
    do_parse!(
        tag!(BAI_MAGIC) >> 
        n_ref: le_i32 >>
        ref_indexes: count!(parse_ref_index, n_ref as usize) >>
        n_no_coor: le_u64 >>
        (
            BaiIndex {
                ref_indexes: ref_indexes,
                n_no_coor: n_no_coor,
            }
        )
    )
);

impl BaiIndex {
    pub fn load<P: AsRef<Path>>(path: P) -> Result<BaiIndex, IndexedReaderError> {
        let mut f = File::open(path)?;

        let mut buf = Vec::new();
        f.read_to_end(&mut buf)?;

        match parse_bai(&buf).to_result() {
            Ok(v) => Ok(v),
            Err(v) => Err(IndexedReaderError::InvalidIndex),
        }
    }
}



/* calculate bin given an alignment covering [beg,end) (zero-based, half-closed-half-open) */
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
/* calculate the list of bins that may overlap with region [beg,end) (zero-based) */
const MAX_BIN: usize = (((1<<18)-1)/7);

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


struct CsiBin {
    bin: u32,
    loffset: u64,
    chunks: Vec<Chunk>,
}

struct CsiRefIndex {
    bins: Vec<CsiBin>,
}

struct CsiIndex {
    min_shift: i32,
    depth: i32,
    aux: Vec<u8>,
    ref_indexes: Vec<CsiRefIndex>,
    n_no_coor: u64,
}

impl CsiIndex {
    pub fn load<P: AsRef<Path>>(path: P) -> Result<CsiIndex, IndexedReaderError> {
        let mut f = File::open(path)?;

        let mut buf = Vec::new();
        f.read_to_end(&mut buf)?;

        match parse_csi(&buf).to_result() {
            Ok(v) => Ok(v),
            Err(v) => Err(IndexedReaderError::InvalidIndex),
        }
    }
}

const CSI_MAGIC: &[u8; 4] = &[b'C', b'S', b'I', 1u8];

named!(parse_csi<&[u8], CsiIndex>,
    do_parse!(
        tag!(CSI_MAGIC) >> 
        min_shift: le_i32 >>
        depth: le_i32 >>
        l_aux: le_i32 >>
        aux: count!(le_u8, l_aux as usize) >>
        n_ref: le_i32 >>
        ref_indexes: count!(parse_csi_ref_index, n_ref as usize) >>
        n_no_coor: le_u64 >>
        (
            CsiIndex {
                min_shift: min_shift,
                depth: depth,
                aux: aux,
                ref_indexes: ref_indexes,
                n_no_coor: n_no_coor,
            }
        )
    )
);


named!(parse_csi_ref_index<&[u8], CsiRefIndex>,
    do_parse!(
        n_bin: le_i32 >>
        bins: count!(parse_csi_bin, n_bin as usize) >>
        (
            CsiRefIndex {
                bins: bins,
            }
        )
    )
);

named!(parse_csi_bin<&[u8], CsiBin>,
    do_parse!(
        bin: le_u32 >>
        loffset: le_u64 >>
        n_chunk: le_i32 >>
        chunks: count!(parse_chunk, n_chunk as usize) >>
        (
            CsiBin {
                bin: bin,
                loffset: loffset,
                chunks: chunks,
            }
        )
    )
);


/* calculate bin given an alignment covering [beg,end) (zero-based, half-close-half-open) */
fn csi_reg2bin(u64 beg, u64 end, i32 min_shift, i32 depth)
{
    let l = depth;
    let s = min_shift,
    let t = ((1<<depth*3) - 1) / 7;

    end -= 1;

    while l > 0 {
        if (beg >> s == end >> s) {
            return t + (beg >> s);
        }

        l -= 1;
        s += 3;
        t -= 1<<l*3;
    }

    return 0;
}


/* calculate the list of bins that may overlap with region [beg,end) (zero-based) */
fn csi_reg2bins(u64 beg, u64 end, i32 min_shift, i32 depth, bins: &mut Vec<u32>) {
    bins.clear();
    let l = 0;
    let n = 0;
    let t = 0;
    let s = min_shift + depth * 3;

    end -= 1;
    while l <= depth {
        let b = t + (beg >> s);
        let e = t + (end >> s);
        bins.extend(b..(e+1))

        s -= 3;
        t += 1 << l*3;
        l += 1;
    }
}

pub trait BgzfIndex {
    fn query(&self, ref_id: usize, beg: u32, end: u32, chunks: &mut Vec<Chunk>);
}

use std::iter::FromIterator;

impl BgzfIndex for BaiIndex {

    /// Fill `chunks` with set of chunks that may contain reads inside the query interval
    fn query(&self, ref_id: usize, beg: u32, end: u32, chunks: &mut Vec<Chunk>) {

        chunks.clear();
        let mut bins = Vec::new();
        reg2bins(beg, end, &mut bins);
        let bins: HashSet<u16> = HashSet::from_iter(bins.into_iter());

        let ref_idx = &self.ref_indexes[ref_id];

        // Linear index is over 16384 (2^14) bp bins
        let linear_index_bin = beg >> (2^14);
        let linear_index_start = ref_idx.ioffsets[linear_index_bin as usize];

        for b in ref_idx.bins.iter() {
            if bins.contains(&(b.bin as u16)) {
                chunks.extend(b.chunks.cloned().filter(|c| c.chunk_end > linear_index_start));
            }
        }

        chunks.sort();
    }
}

impl BgzfIndex for CsiIndex {

    /// Fill `chunks` with set of chunks that may contain reads inside the query interval
    fn query(&self, ref_id: usize, beg: u32, end: u32, chunks: &mut Vec<Chunk>) {

        chunks.clear();
        let mut bins = Vec::new();
        csi_reg2bins(beg, end, self.min_shift, self.depth, &mut bins);
        let bins: HashSet<u16> = HashSet::from_iter(bins.into_iter());

        let ref_idx = &self.ref_indexes[ref_id];

        for b in ref_idx.bins.iter() {
            if bins.contains(&(b.bin as u16)) {
                // FIXME consider linear indexing
                chunks.extend(b.chunks);
            }
        }

        chunks.sort();
    }
}