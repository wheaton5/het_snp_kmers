#[macro_use]
extern crate clap;
extern crate bloom;
extern crate flate2;
extern crate debruijn;
extern crate fnv;

use fnv::FnvHasher;
use std::collections::{HashSet};
use std::hash::BuildHasherDefault;

//type FnvHashMap<K, V> = HashMap<K, V, BuildHasherDefault<FnvHasher>>;
type FnvHashSet<V> = HashSet<V, BuildHasherDefault<FnvHasher>>;

use clap::{App};

use flate2::read::GzDecoder;
use std::io;
use std::io::prelude::*;
use std::fs::File;
use std::io::Write;

use std::cmp::min;

use debruijn::*;
//use debruijn::dna_string::*;
use debruijn::kmer::*;

use std::str;

static mut KMER_SIZE: usize = 0;

use bloom::{ASMS,CountingBloomFilter,BloomFilter};

fn main() {
    let (input_files, min, max, max_error, max_sum, output_hist, estimated_kmers) = load_params();
    let counting_bits: usize = 7;
    let (bloom_kmer_counter, covered_kmers) = count_kmers_fastq(&input_files, counting_bits, estimated_kmers);
    detect_het_kmers(bloom_kmer_counter, &input_files, estimated_kmers, min, max, max_error, max_sum, output_hist.to_string());
}

fn count_kmers_fastq(kmers_in: &Vec<String>, counting_bits: usize, estimated_kmers: u32) -> (CountingBloomFilter, FnvHashSet<u64>) {
    let mut kmer_counts: CountingBloomFilter = CountingBloomFilter::with_rate(counting_bits, 0.05, estimated_kmers);
    let mut coverage_passing_kmers: FnvHashSet<u64> = FnvHashSet::default();
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
                
                //let dna: DnaString = DnaString::from_dna_string(&line.unwrap());
                //for kmer_start in 0..(dna.len() - k_size + 1) {
                for k in KmerX::kmers_from_ascii(&line.unwrap().as_bytes()) {
                    //let k: KmerX = dna.get_kmer(kmer_start);
                    kmer_counts.insert(&min(k.to_u64(), k.rc().to_u64()));
                }
            }
        }
    }
    (kmer_counts, coverage_passing_kmers)
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


fn detect_het_kmers(kmer_counts: CountingBloomFilter, fastqs: &Vec<String>, 
        estimated_kmers: u32, min_coverage: u32, max_coverage: u32, max_error: u32, max_sum: u32, 
        output_hist: String) {
    let mut visited_kmer = BloomFilter::with_rate(0.03, estimated_kmers);
    eprintln!("counting bloom filter created, now second pass to detect het kmers");
    //let max_count: u32 = (2u32).pow(counting_bits as u32) - 1;
    let mut full_hist = match File::create(output_hist) {
            Err(_) => panic!("couldn't open file for writing"),
            Ok(file) => file,
    };
    let middle_base = KX::K()/2;
    
    
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
                //let dna: DnaString = DnaString::from_dna_string(&line.unwrap());
                for mut k in KmerX::kmers_from_ascii(&line.unwrap().as_bytes()) {
                    //for kmer_start in 0..(dna.len() - k_size + 1) {
                    //let mut k: Kmer21 = dna.get_kmer(kmer_start);
                    //let middle_base = k.get(k_size / 2);
                    let krc = k.rc();
                    let to_hash = min(k.to_u64(), krc.to_u64());
                    if !visited_kmer.contains(&to_hash) {
                        visited_kmer.insert(&to_hash);
                        let count1 = kmer_counts.estimate_count(&to_hash);
                        let mut all_alt_count = 0;
                        if count1 > min_coverage && count1 < max_coverage {
                            let alt_bases = get_alts(k.get(middle_base));
                            for (alt_index, alt_base) in alt_bases.iter().enumerate() {
                                k.set_mut(middle_base, *alt_base);
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
                                    k.set_mut(middle_base, *alt_base);
                                    let mut to_write = k.to_string();
                                    to_write.push_str(&"\t");
                                    to_write.push_str(&count1.to_string());
                                    to_write.push_str(&"\t");
                                    k.set_mut(middle_base, *alt_base);
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


type KmerX = VarIntKmer<u64, KX>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct KX;

impl KmerSize for KX {
    fn K() -> usize {
        unsafe {
            KMER_SIZE
        }
    }
}

fn load_params() -> (Vec<String>, u32, u32, u32, u32, String, u32) {
    let yaml = load_yaml!("params.yml");
    let params = App::from_yaml(yaml).get_matches();
    let mut input_files: Vec<String> = Vec::new();
    for input_file in params.values_of("fastqs").unwrap() {
        input_files.push(input_file.to_string());
    }
    let min = params.value_of("min_coverage").unwrap();
    let min: u32 = min.to_string().parse::<u32>().unwrap();
    let max = params.value_of("max_coverage").unwrap();
    let max: u32 = max.to_string().parse::<u32>().unwrap();
    let kmer_size = params.value_of("kmer_size").unwrap_or("21");
    let kmer_size: usize = kmer_size.to_string().parse::<usize>().unwrap();
    unsafe {
        KMER_SIZE = kmer_size;
    }
    let max_error = params.value_of("max_error").unwrap_or("0");
    let max_error: u32 = max_error.to_string().parse::<u32>().unwrap();
    let max_sum = params.value_of("max_total_coverage").unwrap();
    let max_sum: u32 = max_sum.to_string().parse::<u32>().unwrap();
    let output_hist = params.value_of("output_full_hist").unwrap_or("none");
    let estimated_kmers = params.value_of("estimated_kmers").unwrap_or("1000000000");
    let estimated_kmers: u32 = estimated_kmers.to_string().parse::<u32>().unwrap();
    (input_files, min, max, max_error, max_sum, output_hist.to_string(), estimated_kmers)
}
