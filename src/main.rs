#[macro_use]
extern crate clap;
extern crate rayon;
extern crate bloom;
extern crate flate2;
extern crate debruijn;
extern crate fnv;
extern crate dna_io;
use fnv::FnvHasher;
use std::io::Write;
use std::collections::{HashSet};
use std::hash::BuildHasherDefault;

type FnvHashSet<V> = HashSet<V, BuildHasherDefault<FnvHasher>>;
use rayon::prelude::*;

use clap::{App};
use std::sync::Arc;

use std::error::Error;

use std::cmp::min;

use debruijn::*;
use debruijn::kmer::*;

static mut KMER_SIZE: usize = 21;

use bloom::CountingBloomFilter;

fn main() {
    let (input_files, output_hist, parameters) = load_params();
    //let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(parameters.threads).build().unwrap();
    //let mut thread_handles: Vec<std::thread::JoinHandle<Hist>> = Vec::new();

    //for thread_index in 0..parameters.threads {
    //    let t = thread_spawner(input_files.clone(), parameters.clone(), thread_index);
    //    thread_handles.push(t);
    //}

    //wait_and_join_hists(thread_handles, output_hist);
    //let thread_data = thread_pool.install(|| read_and_spawn(input_files, parameters));
    read_and_spawn(input_files, parameters);
}

fn read_and_spawn(input_files: Vec<String>, parameters: Params) {
    let mut all_thread_data = ThreadData::new_set(parameters);
    for input_file in input_files {
        let reader = dna_io::DnaReader::from_path(&input_file);
        for record in reader {
            let seq = Arc::new(record.seq.clone());
            all_thread_data.par_iter_mut().for_each(|thread_data| {
                thread_data.handle_seq(seq.clone());
            });
            //for thread_data in all_thread_data.iter_mut() {
            //    handle_seq(thread_data, record.seq.clone());
            //}
        }
    }
    all_thread_data.par_iter().map(|thread_data| thread_data.detect_het_kmers() ).for_each(drop);
    //for thread_data in all_thread_data.iter() {
    //    thread_data.detect_het_kmers();
    //}
}


fn wait_and_join_hists(thread_handles: Vec<std::thread::JoinHandle<Hist>>, output_hist: String) {
    let mut big_hist: [u32; 1000] = [0; 1000];
    for t in thread_handles {
        let mini_hist = t.join().unwrap();
        for (index, count) in mini_hist.hist.iter().enumerate() {
            big_hist[index] += *count;
        }
    }
    // do something with big_hist
    let path = std::path::Path::new(&output_hist);
    let display = path.display();
    let mut file = match std::fs::File::create(&path) {
        Err(why) => panic!("couldn't create {}: {}",
                           display,
                           why.description()),
        Ok(file) => file,
    };
    for (index, count) in big_hist.iter().enumerate() {
        let to_write = format!("{}\t{}\n",index.to_string(),count.to_string());
        match file.write_all(to_write.as_bytes()) {
            Err(why) => {
                panic!("couldn't write to {}: {}", display,
                                               why.description())
            },
            Ok(_) => (),
        }
    }
}

struct ThreadData {
    kmer_counts: CountingBloomFilter,
    visited_kmers: FnvHashSet<u64>,
    hist: [u32; 1000],
    mod_index: u64,
    total_threads: u64,
    middle_base_mask: u64,
    params: Params,
}

impl ThreadData {
    fn new_set(params: Params) -> Vec<Self> {
        let mut to_ret: Vec<Self> = Vec::new();
        let kmers = params.estimated_kmers / params.threads;
        if kmers > std::u32::MAX as u64 { panic!("cannot support this number of estimated kmers for this number of threads"); }
        for mod_index in 0..params.threads {
            let mut cbf = CountingBloomFilter::with_rate(params.counting_bits, 0.05, kmers as u32);
            let mut visited: FnvHashSet<u64> = FnvHashSet::default();
            let mut hist = [0u32; 1000];
            let mut td = ThreadData{
                kmer_counts: cbf,
                visited_kmers: visited,
                hist: hist,
                mod_index: mod_index,
                total_threads: params.threads,
                middle_base_mask: !( 3 << ( KX::K() - 1 ) ), // make a mask that is 1's outside the two bits at the center of the kmer
                params: params.clone(),
            };
            to_ret.push(td);
        }
        //for ret in &to_ret { println!("mod index {}",ret.mod_index); }
        to_ret
    }

    fn handle_seq(&mut self, seq: Arc<String>) {
        for k in KmerX::kmers_from_ascii(&seq.as_bytes()) {
            let krc = k.rc();
            let middle_base_invariant = min( k.to_u64() & self.middle_base_mask, krc.to_u64() & self.middle_base_mask );
            let mut kmerstring = String::new();
            if middle_base_invariant % self.total_threads == self.mod_index {
                kmerstring.push_str(&KmerX::from_u64(min( k.to_u64(), krc.to_u64() )).to_string());
            }
            println!("thread {}, kmer {} u64 {} mod 8 {} {}", self.mod_index, k.to_string(),k.to_u64(),middle_base_invariant % self.total_threads, kmerstring);
            if middle_base_invariant % self.total_threads != self.mod_index { continue; }
            let to_hash = min( k.to_u64(), krc.to_u64() );
            println!("handle thread {} inserting {}",self.mod_index, KmerX::from_u64(to_hash).to_string());
            match self.kmer_counts.insert_get_count( &to_hash ) {
                a if a >= self.params.min_count => { // add a middle base invariant hash to the coverage_passing_kmers set
                    self.visited_kmers.insert(middle_base_invariant);
                    println!("got above threshold with {} with {}", KmerX::from_u64(to_hash).to_string(), a);
                },
                _ => (),
            }
        }
    }


    fn detect_het_kmers(&self) {
		eprintln!("counting bloom filter created, detecting het kmers");
		//1's except for middle bits representing middle base and upper mask 0'd out
		let middle_base_mask: u64 = !( 3 << ( KX::K() ) ) & KmerX::top_mask( KX::K() );
		let a_mask = 0; // a mask you can | with a middle base invariant kmer to get that kmer with an A in the middle
		let c_mask = ( 1 << ( KX::K() - 1 ) ) | middle_base_mask; // or a C in the middle
		let g_mask = ( 2 << ( KX::K() - 1 ) ) | middle_base_mask; // etc
		let t_mask = ( 3 << ( KX::K() - 1 ) ) | middle_base_mask; // etc
		let masks: [u64; 4] = [ a_mask, c_mask, g_mask, t_mask ];

		let mut full_hist: [u32; 1000] = [0;1000]; // initialize histogram vector
		let mut alt_counts: [u32; 4] = [0,0,0,0]; let mut total_counts = 0;
		let mut best_count = 0; let mut best_count_index = 0;
		let mut second_best_count = 0; let mut second_best_count_index = 0;    
		for middle_base_invariant_kmer in &self.visited_kmers {
			let a_hash = middle_base_invariant_kmer;  // always < than its RC so no check
			let c_hash = middle_base_invariant_kmer | c_mask; // always < its RC so no check
			let g_u64 = middle_base_invariant_kmer | g_mask; //could be > than its RC if palindromic. must check
			let gmer = KmerX::from_u64( g_u64 );
			let g_hash = min( g_u64, gmer.rc().to_u64() );
			let t_u64 = middle_base_invariant_kmer | t_mask; //could be > than its RC if palindromic. must check
			let tmer = KmerX::from_u64( t_u64 );
			let t_hash = min( t_u64, tmer.rc().to_u64() );
            println!("{}",gmer.to_string());
			alt_counts[0] = self.kmer_counts.estimate_count( &a_hash );
			alt_counts[1] = self.kmer_counts.estimate_count( &c_hash );
			alt_counts[2] = self.kmer_counts.estimate_count( &g_hash );
			alt_counts[3] = self.kmer_counts.estimate_count( &t_hash );
            println!("double check {} {} {} {}",alt_counts[0], alt_counts[1], alt_counts[2], alt_counts[3]);
			for (index, count) in alt_counts.iter().enumerate() {
				full_hist[min(full_hist.len()-1,*count as usize)] += 1;
				total_counts += *count;
                println!("\t{}",*count);
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
			if best_count >= self.params.min_count && second_best_count >= self.params.min_count &&
					best_count <= self.params.max_count && total_counts - best_count - second_best_count <= self.params.max_error &&
					best_count + second_best_count <= self.params.max_sum {
				let kmer = KmerX::from_u64(middle_base_invariant_kmer | masks[best_count_index]);
				let kmer2 = KmerX::from_u64(middle_base_invariant_kmer | masks[second_best_count_index]);
				println!("{}\t{}\t{}\t{}",kmer.to_string(), best_count.to_string(), kmer2.to_string(), second_best_count.to_string());
			}
			total_counts = 0; best_count = 0; best_count_index = 0; second_best_count = 0; second_best_count_index = 0; //reset counters
			for i in 0..3 { alt_counts[i] = 0; }
		}
    }
}
fn count_kmers_fastq(kmers_in: Vec<String>, params: Params, thread_index: u64) -> (CountingBloomFilter, FnvHashSet<u64>) {
    let total_kmers = params.estimated_kmers / params.threads;
    if total_kmers > std::u32::MAX as u64 { panic!("cannot deal with this many estimated kmers with this few threads."); }
    let mut kmer_counts: CountingBloomFilter = CountingBloomFilter::with_rate(params.counting_bits, 0.05, total_kmers as u32);
    let mut coverage_passing_kmers: FnvHashSet<u64> = FnvHashSet::default();
    let middle_base_mask: u64 = !( 3 << ( KX::K() - 1 ) ); // make a mask that is 1's outside the two bits at the center of the kmer 
    for kmer_file in &kmers_in {
        let reader = dna_io::DnaReader::from_path(kmer_file);
        for record in reader {
            for k in KmerX::kmers_from_ascii(&record.seq.as_bytes()) {
                let krc = k.rc();
                let middle_base_invariant = min( k.to_u64() & middle_base_mask, krc.to_u64() & middle_base_mask );
                // below is the line that splits thread work and maintains that ...X... and ...Y... kmers end up in the same thread
                if middle_base_invariant % params.threads != thread_index { continue; }
                let to_hash = min( k.to_u64(), krc.to_u64() );
                match kmer_counts.insert_get_count( &to_hash ) {
                    a if a >= params.min_count => { // add a middle base invariant hash to the coverage_passing_kmers set
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

fn detect_het_kmers(kmer_counts: CountingBloomFilter, 
        covered_kmers: FnvHashSet<u64>, params: Params) -> Hist {
    
    eprintln!("counting bloom filter created, detecting het kmers");
    //1's except for middle bits representing middle base and upper mask 0'd out
    let middle_base_mask: u64 = !( 3 << ( KX::K() ) ) & KmerX::top_mask( KX::K() ); 
    let a_mask = 0; // a mask you can | with a middle base invariant kmer to get that kmer with an A in the middle
    let c_mask = ( 1 << ( KX::K() - 1 ) ) | middle_base_mask; // or a C in the middle
    let g_mask = ( 2 << ( KX::K() - 1 ) ) | middle_base_mask; // etc
    let t_mask = ( 3 << ( KX::K() - 1 ) ) | middle_base_mask; // etc
    let masks: [u64; 4] = [ a_mask, c_mask, g_mask, t_mask ];
    
    let _max_count: u32 = (2u32).pow(params.counting_bits as u32) - 1;
    let mut full_hist: [u32; 1000] = [0;1000]; // initialize histogram vector
    let mut alt_counts: [u32; 4] = [0,0,0,0]; let mut total_counts = 0;
    let mut best_count = 0; let mut best_count_index = 0;
    let mut second_best_count = 0; let mut second_best_count_index = 0;
    for middle_base_invariant_kmer in &covered_kmers {
        let a_hash = middle_base_invariant_kmer;  // always < than its RC so no check
        let c_hash = middle_base_invariant_kmer | c_mask; // always < its RC so no check
        let g_u64 = middle_base_invariant_kmer | g_mask; //could be > than its RC if palindromic. must check
        let gmer = KmerX::from_u64( g_u64 );
        let g_hash = min( g_u64, gmer.rc().to_u64() );
        let t_u64 = middle_base_invariant_kmer | t_mask; //could be > than its RC if palindromic. must check
        let tmer = KmerX::from_u64( t_u64 );
        let t_hash = min( t_u64, tmer.rc().to_u64() );
        alt_counts[0] = kmer_counts.estimate_count( &a_hash );
        alt_counts[1] = kmer_counts.estimate_count( &c_hash );
        alt_counts[2] = kmer_counts.estimate_count( &g_hash );
        alt_counts[3] = kmer_counts.estimate_count( &t_hash );
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
        if best_count >= params.min_count && second_best_count >= params.min_count && 
                best_count <= params.max_count && total_counts - best_count - second_best_count <= params.max_error && 
                best_count + second_best_count <= params.max_sum {
            let kmer = KmerX::from_u64(middle_base_invariant_kmer | masks[best_count_index]);
            let kmer2 = KmerX::from_u64(middle_base_invariant_kmer | masks[second_best_count_index]);
            println!("{}\t{}\t{}\t{}",kmer.to_string(), best_count.to_string(), kmer2.to_string(), second_best_count.to_string());
        }
        total_counts = 0; best_count = 0; best_count_index = 0; second_best_count = 0; second_best_count_index = 0; //reset counters
        for i in 0..3 { alt_counts[i] = 0; }
    }
    //let _full_hist = match File::create(output_hist) {
    //        Err(_) => panic!("couldn't open file for writing"),
    //        Ok(file) => file,
    //};
    Hist{ hist: full_hist }
}

struct Hist{
    hist: [u32; 1000],
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

#[derive(Clone)]
struct Params {
    min_count: u32,
    max_count: u32,
    max_error: u32,
    max_sum: u32,
    threads: u64,
    estimated_kmers: u64,
    counting_bits: usize,
}

fn load_params() -> (Vec<String>, String, Params) {
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
    let estimated_kmers: u64 = estimated_kmers.to_string().parse::<u64>().unwrap();
    let threads = params.value_of("threads").unwrap_or("1");
    let threads: u64 = threads.to_string().parse::<u64>().unwrap();
    let counting_bits = params.value_of("counting_bits").unwrap_or("7");
    let counting_bits: usize = counting_bits.to_string().parse::<usize>().unwrap();
    let params: Params = Params{
        min_count: min,
        max_count: max,
        max_error: max_error,
        max_sum: max_sum,
        threads: threads,
        estimated_kmers: estimated_kmers,
        counting_bits: counting_bits,
    };
    (input_files, output_hist.to_string(), params)
}
