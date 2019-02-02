# het_snp_kmers
memory efficient de novo detection of het snp kmers using counting bloom filters.

Install requirements: rust ver 1.3 or later. 
If you do not have rust
```
curl https://sh.rustup.rs -sSf | sh
echo 'export PATH=~/.cargo/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
which cargo
```
For OSX only you will need xz for the htslib dependency. Linux should have the required libraries by default. You can get this with homebrew
```
brew install xz
then find where it installed, and put the include path on your CFLAGS
export CFLAGS=-I/path/to/xz/<version>/include'
or add that to your .bashrc and source it
```
Then you should be able to clone and install the project.
```
git clone git@github.com:wheaton5/het_snp_kmers.git
cd het_snp_kmers
cargo build --release
```

```
./target/release/het_snp_kmers -h
het_snp_kmers 1.0
Haynes Heaton <whheaton@gmail.com>
Finds kmer pairs that are different in the middle base and each have roughly haploid coverage. Meant for illumina data
as an initial step for de novo phasing.

USAGE:
    het_snp_kmers [OPTIONS] --estimated_kmers <estimated_kmers> --inputs <inputs>... --max_coverage <max_coverage> --min_coverage <min_coverage> --output_full_hist <output_full_hist>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
        --estimated_kmers <estimated_kmers>
            estimated total unique kmers. good rule of thumb is roughly 2 * genome size

    -i, --inputs <inputs>...
            input sequence files (fastq,fasta can be gzipped, sam, bam) from which to find het snp kmers

    -k, --kmer_size <kmer_size>                      kmer size to use, defaults to 21
        --max_coverage <max_coverage>                max coverage for each kmer of the pair
        --max_error <max_error>
            max count of the other two kmers with middle base changed. For best results best to be strict and use 0.

        --max_total_coverage <max_total_coverage>    max sum of all kmers with middle base changed
        --min_coverage <min_coverage>                min coverage for each kmer of the pair
        --output_full_hist <output_full_hist>        file name for full kmer histogram
    -t, --threads <threads>                          number of threads to use, defaults to 1
```


Example on made up small test data
```
./target/release/het_snp_kmers --inputs test/data/test.fastq.gz --estimated_kmers 20 --min_coverage 4 --max_coverage 10 --max_total_coverage 80 --max_error 0 --output_full_hist hist.tsv
counting bloom filter created, detecting het kmers
AAAAAAGGGGACCCCCTTTTT   7       AAAAAAGGGGGCCCCCTTTTT   5
```
