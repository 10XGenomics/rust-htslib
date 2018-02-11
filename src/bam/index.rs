use nom::{HexDisplay,Offset,Needed,IResult,ErrorKind,le_u8,le_u16, le_u32, le_i32, le_u64};
use nom::Err;
use std::str::from_utf8;


struct Chunk {
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
        bins: many_m_n!(n_bin as usize, n_bin as usize, parse_bin) >>
        n_intv: le_i32 >>
        ioffsets: many_m_n!(n_intv as usize, n_intv as usize, le_u64) >>
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
        chunks: many_m_n!(n_chunk as usize, n_chunk as usize, parse_chunk) >>
        (
            Bin {
                bin: bin,
                chunks: chunks,
            }
        )
    )
);

const MAGIC: &[u8; 4] = &[b'B', b'A', b'I', 1u8];

named!(parse_bai<&[u8], Index>,
    do_parse!(
        tag!(MAGIC) >> 
        n_ref: le_i32 >>
        ref_indexes: many_m_n!(n_ref as usize, n_ref as usize, parse_ref_index) >>
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
    if (beg>>14 == end>>14) { return ((1<<15)-1)/7 + (beg>>14); }
    if (beg>>17 == end>>17) { return ((1<<12)-1)/7 + (beg>>17); }
    if (beg>>20 == end>>20) { return ((1<<9)-1)/7 + (beg>>20); }
    if (beg>>23 == end>>23) { return ((1<<6)-1)/7 + (beg>>23); }
    if (beg>>26 == end>>26) { return ((1<<3)-1)/7 + (beg>>26); }
    return 0;
}
/* calculate the list of bins that may overlap with region [beg,end) (zero-based) */
const MAX_BIN: u32 = (((1<<18)-1)/7);

pub fn reg2bins(beg: u32, end: u32, list: &mut [u16; MAX_BIN])
{
    let mut i = 0;
    let end = end - 1;

    list[i] = 0;
    i += 1;

    for k in 1 + (beg>>26) .. 1 + (end>>26) {
        list[i] = k;
        i += 1;
    }

    for k in  9 + (beg>>23) ..  9 + (end>>23) {
        list[i] = k;
        i += 1;
    }

    for k in 73 + (beg>>20) .. 73 + (end>>20) {
        list[i] = k;
        i += 1;
    }

    for k in 585 + (beg>>17) .. 585 + (end>>17) {
        list[i] = k;
        i += 1;
    }

    for k in 4681 + (beg>>14) .. 4681 + (end>>14) {
        list[i] = k;
        i += 1;
    }

    return i;
}