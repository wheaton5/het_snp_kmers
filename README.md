# het_snp_kmers
memory efficient de novo detection of het snp kmers using counting bloom filters.

```
./target/release/het_snp_kmers -h
het_snp_kmers
Haynes Heaton <whheaton@gmail.com>
Finds kmer pairs that are different in the middle base and each have roughly haploid coverage. Meant for illumina data
as an initial step for de novo phasing.

USAGE:
    het_snp_kmers [OPTIONS] --estimated_kmers <estimated_kmers> --fastqs <fastqs>... --max_coverage <max_coverage> --max_total_coverage <max_total_coverage> --min_coverage <min_coverage>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
        --estimated_kmers <estimated_kmers>
            estimated total unique kmers. good rule of thumb is roughly 2 * genome size

    -f, --fastqs <fastqs>...                         fastqs from which to find het snp kmers
        --max_coverage <max_coverage>                max coverage for each kmer of the pair
        --max_error <max_error>
            max count of the other two kmers with middle base changed. For best results best to be strict and use 0.

        --max_total_coverage <max_total_coverage>    max sum of all kmers with middle base changed
        --min_coverage <min_coverage>                min coverage for each kmer of the pair
 ```
