use std::error::Error;
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
    start: i32,
    /// stop position of the random barcode
    #[clap(short, long)]
    stop: i32,
}

#[derive(Debug, Deserialize)]
struct BCMapRecord {
    oligo_dt: String,
    rand_hex: String
}

fn parse_bc_map(bc_map: &str) -> HashMap<String, String, ahash::RandomState> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(bc_map).expect("cannot open barcode map file.");

    let s = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
    let mut hm = HashMap::with_hasher(s);

    for result in rdr.deserialize() {
        // Notice that we need to provide a type hint for automatic
        // deserialization.
        let record: BCMapRecord = result.expect("could not deserialize barcode map record.");
        hm.insert(record.oligo_dt, record.rand_hex);
    }
    hm
}

fn main() {
    let opts: Opts = Opts::parse();

    let bcm = parse_bc_map(&opts.bc_map);

    println!("{:?}", bcm);
}
