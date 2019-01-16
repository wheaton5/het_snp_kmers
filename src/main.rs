extern crate clap;
extern crate bloom;
extern crate rust_htslib;
extern crate fnv;
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

use std::env;
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
            .help("fastqs from which to find het snp kmers"))
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
        .get_matches();
    let args: Vec<String> = env::args().collect();
    let input_files: Vec<_> = matches.value_of("fastqs").unwrap().collect();;
    let min = matches.value_of("min_coverage").unwrap();
    let min: u32 = min.to_string().parse::<u32>().unwrap();
    let max = matches.value_of("max_coverage").unwrap();
    let max: u32 = max.to_string().parse::<u32>().unwrap();
    let max_error = matches.value_of("max_error").unwrap_or("0");
    let max_error: u32 = max_error.to_string().parse::<u32>().unwrap();
    let max_sum = matches.value_of("max_total_coverage").unwrap();
    let max_sum: u32 = max_sum.to_string().parse::<u32>().unwrap();
    let estimated_kmers = matches.value_of("estimated_kmers").unwrap_or("1000000000")
    let estimated_kmers: u32 = estimated_kmers.to_string().parse::<u32>().unwrap();
    let k: usize = 21;
    let counting_bits: u32 = 7;
    println!("lets goooooo");    
    let bloom_kmer_counter = count_kmers_fastq(input_filename.to_string(), counting_bits, estimated_kmers, k);
    detect_het_kmers(bloom_kmer_counter, input_filename.to_string(), chrom.to_string(), start, end, k, counting_bits, min, max, max_error, max_sum);
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
                for kmer_start in 0..(dna.len() - k_size) {
                    let k: Kmer21 = dna.get_kmer(kmer_start);
                    kmer_counts.insert_get_count(&min(k.to_u64(), k.rc().to_u64()));
                }
            }
        }
    }
    kmer_counts
}



//fn revcomp(dna: &str) -> String {
//    // result vector
//    let mut rdna: String = String::with_capacity(dna.len()); 
//    // iterate through the input &str
//    for c in dna.chars().rev() {
//        // test the input
//        match is_dna(c) {
//            false => panic!("Input sequence base is not DNA: {}", dna),
//            true => rdna.push(switch_base(c))
//        }
//    }
//    rdna
//}


//fn is_dna(dna: char) -> bool {
//    match dna {
//        'A' | 'a' | 'C' | 'c' | 'G' | 'g' | 'T' | 't' | 'N'| 'n'  => true,
//        _ => false
//    }
//}

//fn switch_base(b: char) -> char {
//    match b {
//        'A' | 'a' => 'T',
//        'C' | 'c' => 'G',
//        'G' | 'g' => 'C',
//        'T' | 't' => 'A',
//        _ => 'N'
//    }
//}

//fn rc_invariant_hash(bytes: &str) -> u64{
//    let rev = revcomp(bytes);
//    let h1 = hash(bytes.as_bytes());
//    let h2 = hash(rev.as_bytes());
//    if h1 < h2 {
//        h1
//    } else {
//        h2
//    }
//}

//fn hash(bytes: &[u8]) -> u64 {
//    let mut hasher = FnvHasher::default();
//    hasher.write(bytes);
//    hasher.finish()
//}

//fn count_kmers_fasta(input_filename: String, chrom: String, start: u32, end: u32, k: usize, counting_bits: usize) -> CountingBloomFilter {
//    let mut kmer_counting_bloom_filter:CountingBloomFilter = CountingBloomFilter::with_rate(counting_bits,0.05,1000000000);
//    
//    let mut reader = Reader::from_path(input_filename).unwrap();
//    while let Some(record) = reader.next() {
//        let record = record.expect("Error reading record");
///        let mut seq: Vec<u8> = record.owned_seq();
//    
//        let bytes &[u8] = &seq;
//        let s = match str::from_utf8(&bytes) {
//            Ok(v) => v,
//            Err(e) => panic!("Invalid UTF-8 sequence: {}", e),
//        };
//        if seq.len() > k {
//            for i in 0..(seq.len()-k+1) {
 //               let h1 = rc_invariant_hash(&s[i..(i+k)]);
                //if "GTGCGCCACCACGCCTGGCCT" == &s[i..(i+k)] || "GTGCGCCACCACGCCTGGCCT" == revcomp(&s[i..(i+k)]) {
                //    println!("searchme {} {}",&s[i..(i+k)], record.pos()+((i+k/2) as i32));
                //}
                //if "CCATCTGTCTCGGCCTCCCAA" == &s[i..(i+k)] || "CCATCTGTCTCGGCCTCCCAA" == revcomp(&s[i..(i+k)]) {
                //    println!("searchme {} {}",&s[i..(i+k)], record.pos()+((i+k/2) as i32));
                //}
//                kmer_counting_bloom_filter.insert_get_count(&h1);
//            }
//        }
//    }
//    kmer_counting_bloom_filter
//}

fn count_kmers_bam(bam_filename: String, chrom: String, start: u32, end: u32, k: usize, counting_bits: usize) -> CountingBloomFilter {
    let mut kmer_counting_bloom_filter:CountingBloomFilter = CountingBloomFilter::with_rate(counting_bits,0.05,1000000000);

    //let mut bam = bam::IndexedReader::from_path(&bam_filename).unwrap();
    //if chrom != "all" {
    //    let tid = bam.header().tid(&chrom.as_bytes()).unwrap();
    //    bam.fetch(tid, start, end).unwrap();
    //}
    let mut bam = bam::Reader::from_path(&bam_filename).unwrap();
    let mut count = 0;
    for r in bam.records() {
        count += 1;
        if count % 100000 == 0 {
            println!("{}",count);
        }
        let record = r.unwrap();
        let seq = record.seq();
        let bytes = seq.as_bytes();
        let s = match str::from_utf8(&bytes) {
            Ok(v) => v,
            Err(e) => panic!("Invalid UTF-8 sequence: {}", e),
        };
        if seq.len() > k {
            for i in 0..(seq.len()-k+1) {
                let h1 = rc_invariant_hash(&s[i..(i+k)]);
                //if "GTGCGCCACCACGCCTGGCCT" == &s[i..(i+k)] || "GTGCGCCACCACGCCTGGCCT" == revcomp(&s[i..(i+k)]) {
                //    println!("searchme {} {}",&s[i..(i+k)], record.pos()+((i+k/2) as i32));
                //}
                //if "CCATCTGTCTCGGCCTCCCAA" == &s[i..(i+k)] || "CCATCTGTCTCGGCCTCCCAA" == revcomp(&s[i..(i+k)]) {
                //    println!("searchme {} {}",&s[i..(i+k)], record.pos()+((i+k/2) as i32));
                //}
                kmer_counting_bloom_filter.insert_get_count(&h1);
            }
        }
    }
    kmer_counting_bloom_filter
}

fn detect_het_kmers(kmer_counts: CountingBloomFilter, bam_filename: String, chrom: String, start: u32, end: u32, 
        k: usize, counting_bits: usize, min: u32, max: u32, max_error: u32, max_sum: u32) {
    let mut visited_kmer = BloomFilter::with_rate(0.03, 1000000000);
    println!("counting bloom filter created, now second pass to detect het kmers");
    let max_count: u32 = (2u32).pow(counting_bits as u32) - 1;
     

    let mut full_hist = match File::create("fullhist.csv") {
        Err(_) => panic!("couldn't open file for writing"),
        Ok(file) => file,
    };
    let mut het_hist = match File::create("hethist.csv") {
        Err(_) => panic!("couldn't open file for writing"),
        Ok(file) => file,
    };
    
    //let mut bam = bam::IndexedReader::from_path(&bam_filename).unwrap();
    //if chrom != "all" {
    //    let tid = bam.header().tid(&chrom.as_bytes()).unwrap();
    //    bam.fetch(tid, start, end).unwrap();
    //}
    let mut bam = bam::Reader::from_path(&bam_filename).unwrap();
    let mut alt_counts =  vec!(0,0,0);
    let mut count = 0;
    for r in bam.records() { 
        count += 1;
        if count % 100000 == 0 {
            println!("{}",count);
        }
        let record = r.unwrap();
        let seq = record.seq();
        let bytes = seq.as_bytes();
        let s = match str::from_utf8(&bytes) {
            Ok(v) => v,
            Err(e) => panic!("Invalid UTF-8 sequence: {}", e),
        };

        if seq.len() > k {
            for i in 0..(seq.len()-k+1) {
                let h1 = rc_invariant_hash(&s[i..(i+k)]);
                if !visited_kmer.contains(&h1) {
                    let count1 = kmer_counts.estimate_count(&h1);
                    let mut all_alt_count = 0;
                    if count1 > min && count1 < max {
                        let alts = get_alts(&s[i..(i+k)]);
                        let mut bailout = false;
                        for (alt_index, alt) in alts.iter().enumerate() {
                            let hash = rc_invariant_hash(&alt);
                            let count2 = kmer_counts.estimate_count(&hash);
                            alt_counts[alt_index] = count2;
                            if count2 == max_count {
                                bailout = true;
                            }
                            if count2 > min {
                                visited_kmer.insert(&hash);
                            }
                            all_alt_count += count2;
                        }
                        if bailout { continue; }

                        for (alt_index, alt) in alts.iter().enumerate() {
                            let count2 = alt_counts[alt_index];                            
                            if count2 > min && count2 < max && all_alt_count - count2 < max_error && count1+count2 <= max_sum {
                                let mut to_write = count2.to_string();
                                to_write.push_str(&"\n");
                                to_write.push_str(&count1.to_string());
                                to_write.push_str(&"\n");
                                het_hist.write_all(to_write.as_bytes()).expect("Unable to write data");
                                println!("{},{},{},{},{},{},{},{},{}",&chrom,
                                    record.pos()+(i as i32)+(k as i32)/2,
                                    &s[i..(i+k)], &alt,count1,count2,alt_counts[0],alt_counts[1],alt_counts[2]);
                                
                            }
                        }
                    }
                    if count1 > 5 {
                        let mut to_write = count1.to_string();
                        to_write.push_str(&"\n");
                        full_hist.write_all(to_write.as_bytes()).expect("Unable to write data");
                    }
                    visited_kmer.insert(&h1);
                }
            }
        }
    }
}


fn get_alts(kmer: &str) -> Vec<String> {
    let mut to_ret: Vec<String> = Vec::new();
    let kmercopy = kmer.as_bytes();
    let mut a = kmer.to_string().into_bytes();
    let mut c = kmer.to_string().into_bytes();
    let mut g = kmer.to_string().into_bytes();
    let mut t = kmer.to_string().into_bytes();
    let index = kmer.len()/2;
    a[index] = b'A';
    c[index] = b'C';
    g[index] = b'G';
    t[index] = b'T';
    let a = str::from_utf8(&a).unwrap().to_string();
    let c = str::from_utf8(&c).unwrap().to_string();
    let g = str::from_utf8(&g).unwrap().to_string();
    let t = str::from_utf8(&t).unwrap().to_string();
    match kmercopy[index] {
        b'A' | b'a' => {
            to_ret.push(c);
            to_ret.push(g);
            to_ret.push(t);
            to_ret
        },
        b'C' | b'c' => {
            to_ret.push(a);
            to_ret.push(g);
            to_ret.push(t);
            to_ret
        },
        b'G' | b'g' => {
            to_ret.push(a);
            to_ret.push(c);
            to_ret.push(t);
            to_ret
        },
        b'T' | b't' => {
            to_ret.push(a);
            to_ret.push(c);
            to_ret.push(g);
            to_ret
        },
        _ => {
            to_ret
        }
    }
}
