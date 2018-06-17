use nom::{le_u8, le_u32, le_i32, le_u64};
use std::path::Path;
use std::fs::File;
use std::io::Read;
use std::iter::FromIterator;

use std::collections::HashSet;
use bam::IndexedReaderError;

use bam::index::reg2bins;

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

        match parse_bai(&buf) {
            Ok((_, value)) => Ok(value),
            Err(_) => Err(IndexedReaderError::InvalidIndex),
        }
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
    min_shift: u32,
    depth: u32,
    aux: Vec<u8>,
    ref_indexes: Vec<CsiRefIndex>,
    n_no_coor: u64,
}

impl CsiIndex {
    pub fn load<P: AsRef<Path>>(path: P) -> Result<CsiIndex, IndexedReaderError> {
        let mut f = File::open(path)?;

        let mut buf = Vec::new();
        f.read_to_end(&mut buf)?;

        match parse_csi(&buf) {
            Ok((_, value)) => Ok(value),
            Err(_) => Err(IndexedReaderError::InvalidIndex),
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
                min_shift: min_shift as u32,
                depth: depth as u32,
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
    fn query(&self, ref_id: usize, beg: u64, end: u64, chunks: &mut Vec<Chunk>);
}

const BAM_MIN_SHIFT: u32 = 3;
const BAM_N_LEVELS: u32 = 5;

impl BgzfIndex for BaiIndex {

    /// Fill `chunks` with set of chunks that may contain reads inside the query interval
    fn query(&self, ref_id: usize, beg: u64, end: u64, chunks: &mut Vec<Chunk>) {

        chunks.clear();
        let mut bins = Vec::new();
        reg2bins(beg, end, BAM_MIN_SHIFT, BAM_N_LEVELS, &mut bins);
        let bins: HashSet<u32> = HashSet::from_iter(bins.into_iter());

        let ref_idx = &self.ref_indexes[ref_id];

        // Linear index is over 16384 (2^14) bp bins
        let linear_index_bin = beg >> (2^14);
        let linear_index_start = ref_idx.ioffsets[linear_index_bin as usize];

        for b in ref_idx.bins.iter() {
            if bins.contains(&b.bin) {
                chunks.extend(b.chunks.iter().cloned().filter(|c| c.chunk_end > linear_index_start));
            }
        }

        chunks.sort();
    }
}

impl BgzfIndex for CsiIndex {

    /// Fill `chunks` with set of chunks that may contain reads inside the query interval
    fn query(&self, ref_id: usize, beg: u64, end: u64, chunks: &mut Vec<Chunk>) {

        chunks.clear();
        let mut bins = Vec::new();
        reg2bins(beg, end, self.min_shift, self.depth, &mut bins);
        let bins: HashSet<u32> = HashSet::from_iter(bins.into_iter());

        let ref_idx = &self.ref_indexes[ref_id];

        for b in ref_idx.bins.iter() {
            if bins.contains(&b.bin) {
                // FIXME consider linear indexing
                chunks.extend(&b.chunks);
            }
        }

        chunks.sort();
    }
}