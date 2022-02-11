# splitp

Streaming SPLiT-seq read pre-processing in Rust.  Currently, this is just a Rust implementation of the Perl pre-processing script by [@jeremymsimon](https://github.com/jeremymsimon) [here](https://github.com/jeremymsimon/SPLITseq).  In the future, capabilities may expand / become more general.

The raison d'etre of the program is simply to perform equivalent processing to the custom Perl script above, but much faster.  To that end, here's a small benchmark on a set of 10,000,000 SPLiT-seq v2 reads.  This is on a 2016 MacBook Pro 2.9 GHz Quad-Core Intel Core i7, 16 GB 2133 MHz LPDDR3 (but both programs are using a single thread).

| program     | runtime     |
| ----------- | ----------- |
| Perl script | 2m 48.6s    |
| splitp      | 5.7s        |
| splitp (pipe to `/dev/null`)      | 3.7s        |

### Usage 

The `splitp` program takes several arguments.  The usage can be printed 
from the command line using `splitp -h`.

```
USAGE:
    splitp [OPTIONS] --read-file <READ_FILE> --bc-map <BC_MAP> --start <START> --end <END>

OPTIONS:
    -b, --bc-map <BC_MAP>          the map of oligo-dT to random hexamers
    -e, --end <END>                end position of the random barcode
    -h, --help                     Print help information
    -o, --one-hamming              consider 1-hamming distance neighbors of random hexamers
    -r, --read-file <READ_FILE>    the input R2 file
    -s, --start <START>            start position of the random barcode
    -V, --version                  Print version information

```

**Please take note** that `splitp` writes the output (processed reads) to stdout, so that the 
output can be directly piped to an input stream of another program (e.g. directly 
to `alevin-fry` via [process substitution](https://tldp.org/LDP/abs/html/process-sub.html)).  This 
also means, if you want to store the processed reads on disk, you can pipe the results directly 
to gzip to compress them e.g.:

```
splitp -r reads.fq -b oligo_hex_bc_mapping.txt -s 87 -e 94 -o | gzip > reads.fq.gz
```

### Notes

The input oligo-dT to random-mer mapping must be provided in a two-column **tab-separated** file.
Further, the first row of the file must be a comment starting with the `#` character.


### Limitations / differences

* Currently, random hexamers can only be searched at a Hamming distance of 0 or 1.
* If the barcode being considered for replacement (BC1) has `N` characters in it, they are all replaced with `A` before lookup
  in the table of random hexamer hamming neighbors.  If there is more than one `N`, this may result in slightly different 
  behavior than the Perl script.
