use std::collections::HashMap;
use std::cmp::max;
use std::io::Write;
use failure::Error;
use byteorder::{LE, WriteBytesExt};

const TBI: u32 = 0;
const CSI: u32 = 1;
const HTS_MIN_MARKER_DIST: u64 = 0x10000;

use bam::index::{hts_bin_bot, hts_bin_parent, hts_bin_first, hts_bin_level, reg2bin};

struct Z {
    last_bin: u32,
    save_bin: u32,
    last_coor: isize,
    last_tid: isize,
    save_tid: isize,
    finished: bool,
    last_off: u64,
    save_off: u64,
    off_beg: u64,
    off_end: u64,
    n_mapped: u64,
    n_unmapped: u64,
}

struct BIdx {
    table: HashMap<u32, Bin>,
}

impl BIdx {
    fn new() -> BIdx {
        BIdx { table: HashMap::new() }
    }

    fn insert(&mut self, bin: u32, beg: u64, end: u64) {
        let entry = self.table.entry(bin).or_insert_with(|| Bin::new());
        entry.pairs.push((beg, end))
    }
}

struct Bin {
    loff: u64,
    pairs: Vec<(u64, u64)>,
}

impl Bin {
    fn new() -> Bin {
        Bin { loff: 0, pairs: Vec::new() }
    }
}

struct LIdx {
    offsets: Vec<u64>,
}

impl LIdx {
    fn new() -> LIdx {
        LIdx {
            offsets: Vec::new(),
        }
    }

    fn insert(&mut self, _beg: isize, _end: isize, offset: u64, min_shift: u32) {
        let beg = (_beg >> min_shift) as usize;
        let end = ((_end - 1) >> min_shift as isize) as usize;
        let n = self.offsets.len();

        for _ in n..beg {
            self.offsets.push(<u64>::max_value());
        }
        for _ in n..=end {
            self.offsets.push(offset);
        }
    }
}

pub struct HtsIdx {
    fmt: u32,
    min_shift: u32,
    n_levels: u32,
    n_bins: u32,
    
    n_no_coor: u64,

    bidx: Vec<BIdx>,
    lidx: Vec<LIdx>,
    meta: Vec<u8>,
    z: Z
}

const EMPTY: isize = 0xffffffffisize;
const EMPU64: u64 = 0xffffffffu64;
const EMPU32: u32 = 0xffffu32;

impl HtsIdx {
    pub fn new(fmt: u32, offset0: u64, min_shift: u32, n_levels: u32) -> HtsIdx {

        let z = Z {
            save_bin: EMPU32,
            save_tid: EMPTY,
            last_tid: EMPTY,
            last_bin: EMPU32,
            save_off: offset0,
            last_off: offset0,
            off_beg:  offset0,
            off_end: offset0,
            last_coor: EMPTY,
            n_mapped: 0,
            n_unmapped: 0,
            finished: false,
        };

        HtsIdx {
            fmt,
            min_shift,
            n_levels,
            n_bins: ((1<<(3 * n_levels + 3)) - 1) / 7,
            n_no_coor: 0,
            bidx: Vec::new(),
            lidx: Vec::new(),
            z,
            meta: Vec::new(),
        }
    }   

    pub fn push(&mut self, tid: isize, mut beg: isize, mut end: isize, offset: u64, is_mapped: bool) {

        let maxpos= 1 << (self.min_shift + self.n_levels * 3);
        if tid < 0 {
            beg = -1;
            end = 0;
        }

        if tid >= 0 && (beg > maxpos || end > maxpos) {
            panic!("position too big: not compatible with index scheme")
        }

        while tid >= self.bidx.len() as isize {
            self.bidx.push(BIdx::new());
            self.lidx.push(LIdx::new());
        }

        if self.z.finished {
            return;
        }

        // change of chromosome
        if self.z.last_tid != tid || (self.z.last_tid >= 0 && tid < 0) {
            if tid >= 0 && self.n_no_coor > 0 {
                panic!("NO_COOR reads are not in a snigle block at the end")
            }
            
            if tid >= 0 && self.bidx[tid as usize].table.len() > 0 {
                panic!("chromosome block not continuous -- data not sorted");
            }

            self.z.last_tid = tid;
            self.z.last_bin = 0xffffffff;
        } else if tid >= 0 && self.z.last_coor > beg {
            panic!("records not position sorted")
        }

        if tid >=0 {
            if is_mapped {
                beg = max(0, beg);
                end = max(1, end);
                // FIXME -- this this last_off or offset?
                self.lidx[tid as usize].insert(beg, end, self.z.last_off, self.min_shift);
            }
        } else {
            self.n_no_coor += 1;
        }

        let bin = reg2bin(beg as u64, end as u64, self.min_shift, self.n_levels);

        if self.z.last_bin != bin {
            let save_tid = self.z.save_tid as usize;
            if self.z.save_bin != EMPU32 {
                self.bidx[save_tid as usize].insert( 
                    self.z.save_bin,
                    self.z.save_off,
                    self.z.last_off)
            }

            // change of chrom
            if self.z.last_bin == EMPU32 && self.z.save_bin == EMPU32 {
                self.z.off_end = self.z.last_off;
                self.bidx[save_tid].insert(
                    self.n_bins + 1,
                    self.z.off_beg,
                    self.z.off_end,
                );

                self.bidx[save_tid].insert(
                    self.n_bins + 1,
                    self.z.n_mapped,
                    self.z.n_unmapped,
                );

                self.z.n_mapped = 0;
                self.z.n_unmapped = 0;
            }

            self.z.save_off = self.z.last_off;
            self.z.save_bin = bin;
            self.z.last_bin = bin;
            self.z.save_tid = tid;
        }

        // Count mapped / unmapped reads
        if is_mapped {
            self.z.n_mapped += 1;
        } else {
            self.z.n_unmapped += 1;
        }

        self.z.last_off = offset;
        self.z.last_coor = beg;
    }

    fn update_loff(&mut self, tid: usize) {

        let bidx = &mut self.bidx[tid];
        let lidx = &mut self.lidx[tid];
        
        // Fill missing values of linear index correctly...
        let offset0 = {
            let metabin = &bidx.table[&(self.n_bins + 1)];
            metabin.pairs[0].0
        };

        let mut i = 0;
        while i < lidx.offsets.len() && lidx.offsets[i] == EMPU64 {
            lidx.offsets[i] = offset0;
            i += 1;
        }

        for (bin, v) in bidx.table.iter_mut() {
            if *bin < self.n_bins {
                let bot_bin = hts_bin_bot(*bin, self.n_levels);
                let new_loff = 
                    if bot_bin < lidx.offsets.len() as u32 {
                        lidx.offsets[bot_bin as usize]
                    } else {
                        0
                    };
                v.loff = new_loff;

            } else {
                // this is a meta-bin
                v.loff = 0;
            }
        } 
    }

    // FIXME
    fn compress_binning(&mut self, tid: usize) {
        let bidx = &mut self.bidx[tid];

        let mut l = self.n_levels;
        while l > 0 {
            let start = hts_bin_first(l);

            for bin in start .. self.n_bins {
                if hts_bin_level(bin) != l { continue }

                let promote_bin = {
                    let v = match bidx.table.get_mut(&bin) { Some(v) => v, _ => continue };
                    v.pairs.sort();

                    // If a bins data all lies within 65kb in the physical file,
                    // push the records up to the parent bin.
                    v.pairs[v.pairs.len() - 1].1 >> 16 - v.pairs[0].0 >> 16 < HTS_MIN_MARKER_DIST 
                };

                if promote_bin {
                    let old_data = bidx.table.remove(&bin).unwrap();
                    let parent_bin = hts_bin_parent(bin);

                    let parent = bidx.table.entry(parent_bin).or_insert_with(|| Bin::new());
                    parent.pairs.extend(old_data.pairs);
                }
            }
            l -= 1;
        }

        // sort the 0 bin
        bidx.table.get_mut(&0).map(|b| b.pairs.sort());

        // Merge adjacent blocks in bins
        for (k, v) in bidx.table.iter_mut() {
            if *k >= self.n_bins { continue; }

            let mut m = 0;
            let list = &mut v.pairs;
            for l in 1 .. list.len() {
                if list[m].1 >> 16 >= list[l].0 >> 16 {
                    if list[m].1 < list[l].1 {
                        list[m].1 = list[l].1
                    }
                } else {
                    m +=1;
                    list[m] = list[l];
                }
            }

            list.truncate(m+1);
        }
    }

    pub fn finish(&mut self, final_offset: u64) {
        // write final metadata blocks
        if self.z.save_tid >= 0 {
            let save_tid = self.z.save_tid as usize;
            let bidx = &mut self.bidx[save_tid];

            bidx.insert(self.z.save_bin, self.z.save_off, final_offset);
            bidx.insert(self.n_bins+1, self.z.off_beg, final_offset);
            bidx.insert(self.n_bins+1, self.z.n_mapped, self.z.n_unmapped);
        }

        for i in 0 .. self.bidx.len() {
            self.update_loff(i);
            self.compress_binning(i);
        }
    }

    pub fn write<W: Write>(&self, w: &mut W) -> Result<(), Error> {
        w.write_i32::<LE>(self.bidx.len() as i32)?;

        if self.fmt == TBI {
            w.write(&self.meta)?;
        }

        for i in 0 .. self.bidx.len() {
            let bidx = &self.bidx[i];
            let lidx = &self.lidx[i];

            w.write_i32::<LE>(bidx.table.len() as i32)?;
            for (k,v) in &bidx.table {
                w.write_u32::<LE>(*k)?;
                if self.fmt == CSI {
                    w.write_u64::<LE>(v.loff)?;
                }

                w.write_i32::<LE>(v.pairs.len() as i32)?;
                for p in &v.pairs {
                    w.write_u64::<LE>(p.0)?;
                    w.write_u64::<LE>(p.1)?;
                }
            }

            if self.fmt != CSI {
                w.write_i32::<LE>(lidx.offsets.len() as i32)?;
                for o in &lidx.offsets {
                    w.write_u64::<LE>(*o)?;
                }
            }
        }

        w.write_u64::<LE>(self.n_no_coor)?;
        Ok(())
    }
}