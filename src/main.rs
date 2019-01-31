#[macro_use]
extern crate clap;
extern crate bloom;
extern crate flate2;
extern crate debruijn;
extern crate fnv;
extern crate dna_io;

use fnv::FnvHasher;
use std::collections::{HashSet};
use std::hash::BuildHasherDefault;

//type FnvHashMap<K, V> = HashMap<K, V, BuildHasherDefault<FnvHasher>>;
type FnvHashSet<V> = HashSet<V, BuildHasherDefault<FnvHasher>>;

use clap::{App};

//use flate2::read::GzDecoder;
//use std::io;
//use std::io::prelude::*;
use std::fs::File;
//use std::io::Write;

use std::cmp::min;

use debruijn::*;
//use debruijn::dna_string::*;
use debruijn::kmer::*;

static mut KMER_SIZE: usize = 0;

use bloom::CountingBloomFilter;

fn main() {
    let (input_files, min, max, max_error, max_sum, output_hist, estimated_kmers) = load_params();
    let counting_bits: usize = 7;
    let (bloom_kmer_counter, covered_kmers) = count_kmers_fastq(&input_files, counting_bits, estimated_kmers, min, max);
    detect_het_kmers(bloom_kmer_counter, covered_kmers, min, max, max_error, max_sum, output_hist.to_string());
}

fn count_kmers_fastq(kmers_in: &Vec<String>, 
                    counting_bits: usize, 
                    estimated_kmers: u32, 
                    min_count: u32, _max_count: u32) -> 
                        (CountingBloomFilter, FnvHashSet<u64>) {
    let mut kmer_counts: CountingBloomFilter = CountingBloomFilter::with_rate(counting_bits, 0.05, estimated_kmers);
    let mut coverage_passing_kmers: FnvHashSet<u64> = FnvHashSet::default();
    let middle_base_mask: u64 = !(3 << (KX::K()-1)); // make a mask that is 1's outside the two bits at the center of the kmer 
    //println!("middle base mask {:#b}",middle_base_mask);
    for kmer_file in kmers_in {
        let reader = dna_io::DnaReader::from_path(kmer_file);
        for record in reader {
            for k in KmerX::kmers_from_ascii(&record.seq.as_bytes()) {
                let to_hash = min(k.to_u64(), k.rc().to_u64());
                match kmer_counts.insert_get_count(&to_hash) {
                    a if a == min_count => { // add a middle base invariant hash to the coverage_passing_kmers set
                        let middle_base_invariant = min(k.to_u64() & middle_base_mask, k.rc().to_u64() & middle_base_mask);
                        coverage_passing_kmers.insert(middle_base_invariant);
                    },
                    //max_count => { // this would save memory but not allow a full histogram to be output which is useful
                    //    let middle_base_invariant = to_hash & middle_base_mask;
                    //    coverage_passing_kmers.remove(&middle_base_invariant);
                    //},
                    _ => (),
                }
            }
        }
    }
    (kmer_counts, coverage_passing_kmers)
}

fn detect_het_kmers(kmer_counts: CountingBloomFilter, covered_kmers: FnvHashSet<u64>, 
        min_count: u32, max_count: u32, max_error: u32, max_sum: u32, 
        output_hist: String) {
    
    eprintln!("counting bloom filter created, detecting het kmers");
    //1's except for middle bits representing middle base and upper mask 0'd out
    let middle_base_mask: u64 = !(3 << (KX::K())) & KmerX::top_mask(KX::K()); 
    let a_mask = 0; // a mask you can | with a middle base invariant kmer to get that kmer with an A in the middle
    let c_mask = (1 << (KX::K()-1)) | middle_base_mask; // or a C in the middle
    let g_mask = (2 << (KX::K()-1)) | middle_base_mask; // etc
    let t_mask = (3 << (KX::K()-1)) | middle_base_mask; // etc
    let masks: [u64; 4] = [a_mask, c_mask, g_mask, t_mask];
    
    //let max_count: u32 = (2u32).pow(counting_bits as u32) - 1;
    let mut full_hist: [u32; 1000] = [0;1000]; // initialize histogram vector
    let mut alt_counts: [u32; 4] = [0,0,0,0]; let mut total_counts = 0;
    let mut best_count = 0; let mut best_count_index = 0;
    let mut second_best_count = 0; let mut second_best_count_index = 0;
    for middle_base_invariant_kmer in &covered_kmers {
        let a_hash = middle_base_invariant_kmer;  // always < than its RC so no check
        let c_hash = middle_base_invariant_kmer | c_mask; // always < its RC so no check
        let g_u64 = middle_base_invariant_kmer | g_mask; //could be > than its RC if palindromic. must check
        let gmer = KmerX::from_u64(g_u64);
        let g_hash = min(g_u64, gmer.rc().to_u64());
        let t_u64 = middle_base_invariant_kmer | t_mask; //could be > than its RC if palindromic. must check
        let tmer = KmerX::from_u64(t_u64);
        let t_hash = min(t_u64, tmer.rc().to_u64());
        alt_counts[0] = kmer_counts.estimate_count(&a_hash);
        alt_counts[1] = kmer_counts.estimate_count(&c_hash);
        alt_counts[2] = kmer_counts.estimate_count(&g_hash);
        alt_counts[3] = kmer_counts.estimate_count(&t_hash);
        for (index, count) in alt_counts.iter().enumerate() {
            full_hist[min(full_hist.len()-1,*count as usize)] += 1;
            total_counts += *count;
            if *count >= best_count {
                second_best_count = best_count;
                second_best_count_index = best_count_index;
                best_count = *count;
                best_count_index = index;
            } else if *count >= second_best_count {
                second_best_count = *count;
                second_best_count_index = index;
            }
        }
        if best_count >= min_count && second_best_count >= min_count && best_count <= max_count && 
            total_counts - best_count - second_best_count <= max_error && best_count + second_best_count <= max_sum {
            let kmer = KmerX::from_u64(middle_base_invariant_kmer | masks[best_count_index]);
            let kmer2 = KmerX::from_u64(middle_base_invariant_kmer | masks[second_best_count_index]);
            println!("{}\t{}\t{}\t{}",kmer.to_string(), best_count.to_string(), kmer2.to_string(), second_best_count.to_string());
        }
        total_counts = 0; best_count = 0; best_count_index = 0; second_best_count = 0; second_best_count_index = 0; //reset counters
        for i in 0..3 { alt_counts[i] = 0; }
    }
    let _full_hist = match File::create(output_hist) {
            Err(_) => panic!("couldn't open file for writing"),
            Ok(file) => file,
    };
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
    for input_file in params.values_of("inputs").unwrap() {
        input_files.push(input_file.to_string());
    }
    let min = params.value_of("min_coverage").unwrap();
    let min: u32 = min.to_string().parse::<u32>().unwrap();
    let max = params.value_of("max_coverage").unwrap();
    let max: u32 = max.to_string().parse::<u32>().unwrap();
    let kmer_size = params.value_of("kmer_size").unwrap_or("21");
    let kmer_size: usize = kmer_size.to_string().parse::<usize>().unwrap();
    if kmer_size % 2 == 0 {
        panic!("kmer size required to be odd");
    }
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
