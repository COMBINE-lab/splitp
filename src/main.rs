use std::error::Error;
use needletail::{parse_fastx_file, Sequence, FastxReader};
use needletail::bitkmer::BitNuclKmer;
use ahash::{AHasher, RandomState};
use csv::{ReaderBuilder, Reader};
use clap::{AppSettings, Parser};
use std::collections::HashMap;
use serde::Deserialize;

// first, reproduce the appproach from
// https://github.com/jeremymsimon/SPLITseq/blob/main/Preprocess_SPLITseq_collapse_bcSharing.pl


/// This doc string acts as a help message when the user runs '--help'
/// as do all doc strings on fields
#[derive(Parser)]
#[clap(version = "1.0", author = "Rob P. <rob@cs.umd.edu>")]
struct Opts {
    /// the input R2 file
    #[clap(short, long)]
    read_file: String,
    /// the map of oligo-dT to random hexamers
    #[clap(short, long)]
    bc_map: String,
    /// start position of the random barcode
    #[clap(short, long)]
    start: usize,
    /// stop position of the random barcode
    #[clap(short, long)]
    stop: usize,
}

#[derive(Debug, Deserialize)]
struct BCMapRecord {
    oligo_dt: String,
    rand_hex: String
}

fn parse_bc_map(bc_map: &str) -> HashMap<u64, u64, ahash::RandomState> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .comment(Some(b'#'))
        .has_headers(false)
        .from_path(bc_map).expect("cannot open barcode file");

    let s = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
    let mut hm = HashMap::with_hasher(s);

    for result in rdr.deserialize() {
        // Notice that we need to provide a type hint for automatic
        // deserialization.
        let record: BCMapRecord = result.expect("could not deserialize barcode map record.");
        
        if let Some((_, rh, _)) = BitNuclKmer::new(record.rand_hex.as_bytes(), record.rand_hex.len() as u8, false).next() {
            if let Some((_, odt, _)) = BitNuclKmer::new(record.oligo_dt.as_bytes(), record.oligo_dt.len() as u8, false).next() {
                hm.insert(rh.0, odt.0);
            }
        }
    }
    hm
}

fn main() {
    let opts: Opts = Opts::parse();

    let bcm = parse_bc_map(&opts.bc_map);
    let bclen = (opts.stop - opts.start + 1) as u8;

    let mut reader = parse_fastx_file(&opts.read_file).expect("valid path/file");    
    while let Some(record) = reader.next() {
        let seqrec = record.expect("invalid record");
        let mut seq = seqrec.seq().into_owned();
        let fivep = &seq[..opts.start];
        let threep = &seq[opts.stop..];

        if let Some((_, bc, _)) = BitNuclKmer::new(&seq[opts.start..opts.start+(bclen as usize)], bclen, false).next() {
            //println!("{:?}", String::from_utf8(seq[opts.start..opts.stop].to_vec()));
            if let Some(oligo_dt) = bcm.get(&bc.0) {
                println!("found it");
            }
        }
    }

    println!("{:?}", bcm);
}
