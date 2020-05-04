#[macro_use]
extern crate clap;
extern crate fnv;
extern crate hashbrown;
extern crate rand;
extern crate debruijn;
use debruijn::*;
use debruijn::kmer::*;

use rand::Rng;
use rand::rngs::StdRng;
use rand::SeedableRng;

use std::io::BufReader;
use std::io::BufRead;
use std::fs::File;
use std::io::Write;
use std::error::Error;
use std::str;

use hashbrown::{HashMap,HashSet};
use fnv::FnvHasher;
use std::hash::BuildHasherDefault;
type FnvHashMap<T,V> = HashMap<T,V, BuildHasherDefault<FnvHasher>>;
type FnvHashSet<T> = HashSet<T, BuildHasherDefault<FnvHasher>>;

use clap::{App};

use std::cmp::min;

fn main() {
    let params = load_params();
    let (het_kmers, het_kmer_pairs) = load_kmers(&params);
}


fn load_kmers(params: &Params) -> (FnvHashSet<Vec<u8>>, FnvHashMap<Vec<u8>, Vec<u8>>) {
    let mut kmer_counts: FnvHashMap<Vec<u8>, [u16; 4]> = FnvHashMap::default();
    let mut set_to_ret: FnvHashSet<Vec<u8>> = FnvHashSet::default();
    let mut map_to_ret: FnvHashMap<Vec<u8>, Vec<u8>> = FnvHashMap::default();
    let seed: [u8; 32] = [4; 32]; // guaranteed random number determined by fair dice roll https://xkcd.com/221/
    let mut rng: StdRng = SeedableRng::from_seed(seed);
    
    let f = File::open(params.kmer_counts_file.to_string()).expect("Unable to open kmer counts file");
    let mut f = BufReader::new(f);
    let mut base_before_middle = b'0';
    let mut index_before_middle: usize = 0;
    let mut middle_index: usize = 0;
    let mut buf = vec![];
    let mut buf2 = vec![];
    let mut index = 0;
    let mut total_counts = 0; let mut best_count = 0; let mut best_count_index = 0;
    let mut second_best_count = 0; let mut second_best_count_index = 0;
    let base_index = [b'A', b'C', b'G', b'T'];
    loop {
        let num_bytes = f.read_until(b'\t', &mut buf).expect("failure to read kmer");
        if num_bytes == 0 { break; }
        if index == 0 {
            middle_index = (num_bytes - 1)/2;
            index_before_middle = middle_index - 1;
        }
        if buf[index_before_middle] != base_before_middle {
            // find kmer pairs in this set
            for (invariant, counts) in &kmer_counts {
                for (index, count) in counts.iter().enumerate() {
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
                //println!("{} {} {} {}", total_counts, best_count, second_best_count, params.max_error);
                if best_count >= params.min_count && second_best_count >= params.min_count && //PAIRED HET
                        best_count <= params.max_count && total_counts - best_count - second_best_count <= params.max_error &&
                        best_count + second_best_count <= params.max_sum {
                    let mut kmer = invariant.clone();
                    kmer[middle_index] = base_index[best_count_index];
                    let mut kmer2 = invariant.clone();
                    kmer2[middle_index] = base_index[second_best_count_index];
                    if rng.gen_range(0,2) == 0 {
                        println!("{}\t{}\t{}\t{}",str::from_utf8(&kmer).unwrap(), 
                            best_count, str::from_utf8(&kmer2).unwrap(), second_best_count);
                    } else {
                        println!("{}\t{}\t{}\t{}",str::from_utf8(&kmer2).unwrap(), 
                            second_best_count, str::from_utf8(&kmer).unwrap(), second_best_count);
                    }
                    set_to_ret.insert(kmer.clone());
                    //set_to_ret.insert(kmer2.clone());
                    map_to_ret.insert(kmer.clone(), kmer2.clone());
                    map_to_ret.insert(kmer2, kmer);
                } else if best_count >= params.max_count && best_count <= params.max_sum && 
                    total_counts - best_count <= params.max_error { // HOMOZYGOUS
                    let mut kmer = invariant.clone();//KmerX::from_u64(middle_base_invariant_kmer | masks[best_count_index]);
                    kmer[middle_index] = base_index[best_count_index];
                    let kmer_u64 = twobit(&kmer);
                    if kmer_u64 % params.hom_modimizer == 0 {
                        println!("{}\t{}\tHOM\t.",str::from_utf8(&kmer).unwrap(), best_count);
                    }
                } else if best_count >= params.min_count && best_count <= params.max_count &&
                    total_counts - best_count <= params.max_error { // UNPAIRED HET
                    
                    let mut kmer = invariant.clone();//KmerX::from_u64(middle_base_invariant_kmer | masks[best_count_index]);
                    kmer[middle_index] = base_index[best_count_index];
                    let kmer_u64 = twobit(&kmer);
                    //if kmer_u64 % params.hom_modimizer == 0 {
                    //    println!("{}\t{}\tHET\t.",str::from_utf8(&kmer).unwrap(), best_count);
                    //}
                }
                total_counts = 0; best_count = 0; best_count_index = 0; second_best_count = 0; second_best_count_index = 0; //reset counters
            }
            kmer_counts.drain();
            base_before_middle = buf[index_before_middle];
        }
        
        let num_bytes2 = f.read_until(b'\n', &mut buf2).expect("failure to read count");
        if num_bytes2 == 0 { break; }
        let count: u16 = str::from_utf8(&buf2[0..(num_bytes2-1)]).unwrap().parse::<u16>().unwrap();
        let middle_base = buf[index_before_middle + 1];
        let mut middle_base_invariant = buf[0..(num_bytes-1)].to_vec();
        middle_base_invariant[index_before_middle + 1] = b'N';
        let counts = kmer_counts.entry(middle_base_invariant).or_insert([0; 4]);
        let countdex = match middle_base {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => 0,
        };
        counts[countdex] = count;       
        buf.clear();
        buf2.clear();
        index += 1;
    }
    (set_to_ret, map_to_ret)
}

fn twobit(kmer: &[u8]) -> u64 {
    let mut kmer_u64: u64 = 0;
    for base in kmer {
        kmer_u64 = kmer_u64 << 2;
        let base: u64 = match base {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => { assert!(false); 0},
        };
        //println!("base {}",base);
        assert!(base < 4);
        kmer_u64 = kmer_u64 | (base as u64);
    }
    kmer_u64
}

#[derive(Clone)]
struct Params {
    kmer_counts_file: String, 
    min_count: u16,
    max_count: u16,
    max_error: u16,
    max_sum: u16,
    counting_bits: usize,
    threads: usize,
    estimated_kmers: u64,
    unpaired_het_modimizer: u64,
    hom_modimizer: u64,
}

fn load_params() -> Params {
    let yaml = load_yaml!("params.yml");
    let params = App::from_yaml(yaml).get_matches();
    let kmer_counts = params.value_of("kmer_counts").unwrap();
    let min = params.value_of("min_coverage").unwrap();
    let min: u16 = min.to_string().parse::<u16>().unwrap();
    let max = params.value_of("max_coverage").unwrap();
    let max: u16 = max.to_string().parse::<u16>().unwrap();
    let max_error = params.value_of("max_error").unwrap_or("0");
    let max_error: u16 = max_error.to_string().parse::<u16>().unwrap();
    let max_sum = params.value_of("max_total_coverage").unwrap();
    let max_sum: u16 = max_sum.to_string().parse::<u16>().unwrap();
    let unpaired_het_modimizer = params.value_of("unpaired_het_modimizer").unwrap();
    let unpaired_het_modimizer: u64 = unpaired_het_modimizer.to_string().parse::<u64>().unwrap();
    let hom_modimizer = params.value_of("hom_modimizer").unwrap();
    let hom_modimizer: u64 = hom_modimizer.to_string().parse::<u64>().unwrap();
    Params{
        kmer_counts_file: kmer_counts.to_string(),
        min_count: min,
        max_count: max,
        max_error: max_error,
        max_sum: max_sum,
        counting_bits: 0,
        estimated_kmers: 0,
        threads: 0,
        unpaired_het_modimizer: unpaired_het_modimizer,
        hom_modimizer: hom_modimizer,
    }
}
