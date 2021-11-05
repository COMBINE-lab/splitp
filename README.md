# splitp

Streaming SPLiT-seq read pre-processing in Rust.  Currently, this is just a Rust implementation of the Perl pre-processing script by [@jeremymsimon](https://github.com/jeremymsimon) [here](https://github.com/jeremymsimon/SPLITseq).  In the future, capabilities may expand / become more general.

The raison d'etre of the program is simply to perform equivalent processing to the custom Perl script above, but much faster.  To that end, here's a small benchmark on a set of 10,000,000 SPLiT-seq v2 reads.  This is on a 2016 MacBook Pro 2.9 GHz Quad-Core Intel Core i7, 16 GB 2133 MHz LPDDR3 (but both programs are using a single thread).

| program     | runtime     |
| ----------- | ----------- |
| Perl script | 2m 48.6s    |
| splitp      | 5.7s        |
| splitp (pipe to `/dev/null`)      | 3.7s        |

### Limitations / differences

* Currently, random hexamers can only be searched at a Hamming distance of 0 or 1.
* If the barcode being considered for replacement (BC1) has `N` characters in it, they are all replaced with `A` before lookup
  in the table of random hexamer hamming neighbors.  If there is more than one `N`, this may result in slightly different 
  behavior than the Perl script.
