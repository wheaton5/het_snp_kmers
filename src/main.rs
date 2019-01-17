extern crate clap;
extern crate bloom;
extern crate flate2;
extern crate debruijn;

use clap::{Arg, App};

use flate2::read::GzDecoder;
use std::io;
use std::io::prelude::*;
use std::fs::File;
use std::io::Write;

use std::cmp::min;

use debruijn::*;
use debruijn::dna_string::*;
use debruijn::kmer::*;

use std::str;

use bloom::{ASMS,CountingBloomFilter,BloomFilter};

fn main() {
    let matches = App::new("het_snp_kmers")
        .author("Haynes Heaton <whheaton@gmail.com>")
        .about("Finds kmer pairs that are different in the middle base and each have roughly haploid coverage. Meant for illumina data as an initial step for de novo phasing.")
        .arg(Arg::with_name("fastqs")
            .long("fastqs")
            .short("f")
            .takes_value(true)
            .multiple(true)
            .required(true)
            .help("gzipped fastqs from which to find het snp kmers"))
        .arg(Arg::with_name("min_coverage")
            .long("min_coverage")
            .takes_value(true)
            .required(true)
            .help("min coverage for each kmer of the pair"))
        .arg(Arg::with_name("max_coverage")
            .long("max_coverage")
            .takes_value(true)
            .required(true)
            .help("max coverage for each kmer of the pair"))
        .arg(Arg::with_name("max_error")
            .long("max_error")
            .takes_value(true)
            .help("max count of the other two kmers with middle base changed. For best results best to be strict and use 0."))
        .arg(Arg::with_name("max_total_coverage")
            .long("max_total_coverage")
            .takes_value(true)
            .required(true)
            .help("max sum of all kmers with middle base changed"))
        .arg(Arg::with_name("estimated_kmers")
            .long("estimated_kmers")
            .required(true)
            .takes_value(true)
            .help("estimated total unique kmers. good rule of thumb is roughly 2 * genome size"))
        .arg(Arg::with_name("output_full_hist")
            .long("output_full_hist")
            .required(true)
            .takes_value(true)
            .help("output full kmer count histogram to this file if specified otherwise does not output histogram"))
        .get_matches();
    let input_files: Vec<_> = matches.values_of("fastqs").unwrap().collect();;
    let min = matches.value_of("min_coverage").unwrap();
    let min: u32 = min.to_string().parse::<u32>().unwrap();
    let max = matches.value_of("max_coverage").unwrap();
    let max: u32 = max.to_string().parse::<u32>().unwrap();
    let max_error = matches.value_of("max_error").unwrap_or("0");
    let max_error: u32 = max_error.to_string().parse::<u32>().unwrap();
    let max_sum = matches.value_of("max_total_coverage").unwrap();
    let max_sum: u32 = max_sum.to_string().parse::<u32>().unwrap();
    let output_hist = matches.value_of("output_full_hist").unwrap_or("none");
    let estimated_kmers = matches.value_of("estimated_kmers").unwrap_or("1000000000");
    let estimated_kmers: u32 = estimated_kmers.to_string().parse::<u32>().unwrap();
    let k: usize = 21;
    let counting_bits: usize = 7;
    let bloom_kmer_counter = count_kmers_fastq(&input_files, counting_bits, estimated_kmers, k);
    detect_het_kmers(bloom_kmer_counter, &input_files, k, estimated_kmers, min, max, max_error, max_sum, output_hist.to_string());
}

fn count_kmers_fastq(kmers_in: &Vec<&str>, counting_bits: usize, estimated_kmers: u32, k_size: usize) -> CountingBloomFilter {
    let mut kmer_counts: CountingBloomFilter = CountingBloomFilter::with_rate(counting_bits, 0.05, estimated_kmers);
    for kmer_file in kmers_in {
        let file = match File::open(kmer_file) {
            Ok(file) => file,
            Err(error) => {
                panic!("There was a problem opening the file: {:?}", error)
            },
        };
        let gz = GzDecoder::new(file);
        for (line_number, line) in io::BufReader::new(gz).lines().enumerate() {
            
            if line_number % 4 == 1 {
                
                let dna: DnaString = DnaString::from_dna_string(&line.unwrap());
                for kmer_start in 0..(dna.len() - k_size + 1) {
                    let k: Kmer21 = dna.get_kmer(kmer_start);
                    kmer_counts.insert(&min(k.to_u64(), k.rc().to_u64()));
                }
            }
        }
    }
    kmer_counts
}

const A_ALT: [u8; 3] = [1,2,3];
const C_ALT: [u8; 3] = [0,2,3];
const G_ALT: [u8; 3] = [0,1,3];
const T_ALT: [u8; 3] = [0,1,2];

fn get_alts(current_base: u8) -> [u8; 3] {
    match current_base {
        0 => A_ALT,
        1 => C_ALT,
        2 => G_ALT,
        3 => T_ALT,
        _ => panic!("not dna"),
    }
} 


fn detect_het_kmers(kmer_counts: CountingBloomFilter, fastqs: &Vec<&str>, 
        k_size: usize, estimated_kmers: u32, min_coverage: u32, max_coverage: u32, max_error: u32, max_sum: u32, 
        output_hist: String) {
    let mut visited_kmer = BloomFilter::with_rate(0.03, estimated_kmers);
    eprintln!("counting bloom filter created, now second pass to detect het kmers");
    //let max_count: u32 = (2u32).pow(counting_bits as u32) - 1;
    let mut full_hist = match File::create(output_hist) {
            Err(_) => panic!("couldn't open file for writing"),
            Ok(file) => file,
    };
    
    
    let mut alt_counts =  vec!(0,0,0);
    //let mut count = 0;
    for kmer_file in fastqs {
        let file = match File::open(kmer_file) {
            Ok(file) => file,
            Err(error) => {
                panic!("There was a problem opening file {} with error {}", kmer_file, error)
            },
        };
        let gz = GzDecoder::new(file);
        for (line_number, line) in io::BufReader::new(gz).lines().enumerate() {
            if line_number % 4 == 1 {
                let dna: DnaString = DnaString::from_dna_string(&line.unwrap());
                for kmer_start in 0..(dna.len() - k_size + 1) {
                    let mut k: Kmer21 = dna.get_kmer(kmer_start);
                    let middle_base = k.get(k_size / 2);
                    let krc = k.rc();
                    let to_hash = min(k.to_u64(), krc.to_u64());
                    if !visited_kmer.contains(&to_hash) {
                        visited_kmer.insert(&to_hash);
                        let count1 = kmer_counts.estimate_count(&to_hash);
                        let mut all_alt_count = 0;
                        if count1 > min_coverage && count1 < max_coverage {
                            let alt_bases = get_alts(k.get(k_size / 2));
                            for (alt_index, alt_base) in alt_bases.iter().enumerate() {
                                k.set_mut(k_size / 2, *alt_base);
                                let krc = k.rc();
                                let to_hash = min(k.to_u64(), krc.to_u64());
                                let count2 = kmer_counts.estimate_count(&to_hash);
                                alt_counts[alt_index] = count2;
                                if count2 > min_coverage {
                                    visited_kmer.insert(&to_hash);
                                }
                                all_alt_count += count2;
                            }
                            for (alt_index, alt_base) in alt_bases.iter().enumerate() {
                                let count2 = alt_counts[alt_index];
                                //println!("{} > {} && {} < {} && {} < {} && {} <= {}",count2,min_coverage, count2, max_coverage,all_alt_count-count2, max_error, count1+count2, max_sum);
                                if count2 > min_coverage && count2 < max_coverage && all_alt_count - count2 <= max_error && count1 + count2 <= max_sum {
                                    //println!(" got here ");
                                    k.set_mut(k_size/2, middle_base);
                                    let mut to_write = k.to_string();
                                    to_write.push_str(&"\t");
                                    to_write.push_str(&count1.to_string());
                                    to_write.push_str(&"\t");
                                    k.set_mut(k_size/2, *alt_base);
                                    to_write.push_str(&k.to_string());
                                    to_write.push_str(&"\t");
                                    to_write.push_str(&count2.to_string());
                                    to_write.push_str(&"\n");
                                    //het_hist.write_all(to_write.as_bytes()).expect("Unable to write data");
                                    //ddif let Some(ref mut x) = &het_hist {
                                    //    x.write_all(to_write.as_bytes()).expect("Unable to write data");
                                    //} else {
                                    println!("{}",to_write);
                                    //}
                                    //match het_hist {
                                    //    Some(mut x) => x.write_all(to_write.as_bytes()).expect("Unable to write data"),
                                    //    None => println!("{}",to_write),
                                    //}
                                }
                            }
                        }
                        if count1 > 5 {
                            let mut to_write = count1.to_string();
                            to_write.push_str(&"\n");
                            //full_hist.write_all(to_write.as_bytes()).expect("Unable to write data");
                            full_hist.write_all(to_write.as_bytes()).expect("Unable to write data");
                        }
                    }
                }
            }
        }
    }
}


type Kmer21 = VarIntKmer<u64, K21>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct K21;

impl KmerSize for K21 {
    fn K() -> usize {
        21
    }
}
