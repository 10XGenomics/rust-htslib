use nom::{le_u8, le_u32, le_i32, le_u64};

use std::collections::HashSet;

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

struct Index {
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

named!(parse_bai<&[u8], Index>,
    do_parse!(
        tag!(BAI_MAGIC) >> 
        n_ref: le_i32 >>
        ref_indexes: count!(parse_ref_index, n_ref as usize) >>
        n_no_coor: le_u64 >>
        (
            Index {
                ref_indexes: ref_indexes,
                n_no_coor: n_no_coor,
            }
        )
    )
);


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

pub trait BgzfIndex {
    fn chunks(&self, ref_id: usize, beg: u32, end: u32, chunks: &mut Vec<Chunk>);
}

use std::iter::FromIterator;

impl BgzfIndex for Index {
    fn chunks(&self, ref_id: usize, beg: u32, end: u32, chunks: &mut Vec<Chunk>) {
        let mut bins = Vec::new();
        reg2bins(beg, end, &mut bins);
        let bins: HashSet<u16> = HashSet::from_iter(bins.into_iter());

        let ref_idx = &self.ref_indexes[ref_id];

        for b in ref_idx.bins.iter() {
            if bins.contains(&(b.bin as u16)) {

                
            }
        }


    }
}