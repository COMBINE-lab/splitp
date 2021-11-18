use clap::Parser;
use needletail::bitkmer::BitNuclKmer;
use needletail::parse_fastx_file;
use serde::Deserialize;
use std::collections::HashMap;
use std::io::Write;
use std::str;
mod utils;

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
    /// end position of the random barcode
    #[clap(short, long)]
    end: usize,
    /// consider 1-hamming distance neighbors of random hexamers
    #[clap(short, long)]
    one_hamming: bool,
}

#[derive(Debug, Deserialize)]
struct BCMapRecord {
    oligo_dt: String,
    rand_hex: String,
}


/// Parse a 2-column tab-separated file (given by `bc_map`) that contains the 
/// mapping of oligo-dT sequences to random-mer sequences.  The file 
/// *must* be tab-separated, and *must* begin with a header line starting 
/// with `#`.
///
/// If the `one_edit` argument is `true` then all one-edit neighbors of an random-mer will 
/// be mapped to the corresponding oligo-dT.  Otherwise, only exactly matching random-mers 
/// will be matched to the corresponding oligo-dT.
fn parse_bc_map(bc_map: &str, one_edit: bool) -> HashMap<u64, String, ahash::RandomState> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .comment(Some(b'#'))
        .has_headers(false)
        .from_path(bc_map)
        .expect("cannot open barcode file");

    let s = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
    let mut hm = HashMap::with_hasher(s);
    let s2 = ahash::RandomState::with_seeds(2u64, 7u64, 1u64, 8u64);
    let mut hs = HashMap::with_hasher(s2);
    let mut collisions = 0usize;
    let mut neighbor_vec = Vec::new();

    for result in rdr.deserialize() {
        // Notice that we need to provide a type hint for automatic
        // deserialization.
        let record: BCMapRecord = result.expect("could not deserialize barcode map record.");

        if let Some((_, rh, _)) = BitNuclKmer::new(
            record.rand_hex.as_bytes(),
            record.rand_hex.len() as u8,
            false,
        )
        .next()
        {
            let count = hs.entry(rh.0).or_insert(0usize);
            if *count > 0{
                collisions += 1;
            }
            *count += 1;
       
            neighbor_vec.push((rh.0, record.oligo_dt.clone()));
            // if we are considering one-edit neighbors, then generate
            // of the neighbors of the random-mer, and prepare them 
            // for insertion into the map with the corresponding 
            // oligo-dT.
            if one_edit {
                for n in utils::get_all_snps(rh.0, record.rand_hex.len()) {
                    let count = hs.entry(n).or_insert(0usize);
                    if *count > 0 {
                        collisions += 1;
                    }
                    *count += 1;
                    neighbor_vec.push((n, record.oligo_dt.clone()));
                }
            }
        }
    }
    // how many collisions do we observe between the 1-edit
    // neighbors of the random hexamers?
    if collisions > 0 {
        eprintln!(
            "Observed {} collisions in the edit distance neighbors of the random hexamers.",
            collisions
        );
    }
    for (k, v) in neighbor_vec {
        if let Some(count) = hs.get(&k) {     
            if *count == 1 {
                hm.insert(k, v);
            } else {
                eprintln!("{} different random mers were within 1-hamming of {}", count, k);
            }
        }
    }
    hm
}

fn main() {
    // parse the options
    let opts: Opts = Opts::parse();

    // read in the random hexamer to oligo-dT map
    // if edit distance bound is 0, only populate with
    // exact matches.  Otherwise, populate with all
    // 1-edit neighbors of random hexamers.
    let bcm = parse_bc_map(&opts.bc_map, opts.one_hamming);

    // the barcode length
    let bclen = ((opts.end - opts.start) + 1) as u8;

    // lock stdout and buffer so we can write to it quickly.
    let stdout = std::io::stdout();
    let lock = stdout.lock();
    let mut buf = std::io::BufWriter::with_capacity(32 * 1024, lock);

    // for now, we're assuming FASTQ and not FASTA.
    let mut reader = parse_fastx_file(&opts.read_file).expect("valid path/file");

    while let Some(record) = reader.next() {
        let seqrec = record.expect("invalid record");
        let seq = seqrec.seq().into_owned();
        // get the parts of the read:
        // fivep = before the barcode
        // bcs = the barcode
        // threep = after the barcode
        let fivep = &seq[..opts.start - 1];
        let bcs = &seq[opts.start - 1..opts.start - 1 + (bclen as usize)];
        let threep = &seq[opts.end..];

        // the fallthrough numeric encoding of
        // the barcode
        let mut bc = 0;

        // currently, this is the only tricky bit.
        // try to convert the barcode to a bitkmer.
        // this will fail only if there are non-ACGT nucleotides in it.
        // in that case, try and replace any `N` with an `A` and try again.
        // if it **still** fails for some reason, just fall back to the
        // fake 0 barcode encoding.
        match BitNuclKmer::new(bcs, bclen, false).next() {
            Some((_, bck, _)) => {
                bc = bck.0;
            }
            None => {
                let sbcs = unsafe { str::from_utf8_unchecked(bcs) };
                let nbcs = str::replace(sbcs, "N", "A");
                if let Some((_, bck, _)) = BitNuclKmer::new(nbcs.as_bytes(), bclen, false).next() {
                    bc = bck.0;
                }
            }
        }

        // write the read info for the part before the
        // sequence.
        buf.write_all(b"@").unwrap();
        buf.write_all(seqrec.id()).unwrap();
        buf.write_all(b"\n").unwrap();

        // try to lookup the barcode in the translation & correction
        // map.  If we find a hit, then write the corrected read
        // otherwise, wrtie the original read.
        if let Some(oligo_dt) = bcm.get(&bc) {
            /*eprintln!(
                "{} => {}",
                String::from_utf8(seq[opts.start - 1..opts.start - 1 + (bclen as usize)].to_vec())
                    .unwrap(),
                oligo_dt
            );
            */
            buf.write_all(fivep).unwrap();
            buf.write_all(oligo_dt.as_bytes()).unwrap();
            buf.write_all(threep).unwrap();
        } else {
            buf.write_all(&seqrec.seq()).unwrap();
        }
        buf.write_all(b"\n+").unwrap();
        buf.write_all(seqrec.id()).unwrap();
        buf.write_all(b"\n").unwrap();
        buf.write_all(seqrec.qual().unwrap()).unwrap();
        buf.write_all(b"\n").unwrap();
    }
    buf.flush().unwrap();
}
