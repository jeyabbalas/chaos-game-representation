use std::collections::HashMap;
use std::convert::TryInto;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::Path;
use flate2::read::GzDecoder;
use rand::prelude::ThreadRng;
use rand::{
    distributions::{Distribution, Standard, WeightedIndex},
    Rng,
    thread_rng, 
};


pub fn convert_decimal_to_binary(decimal: u32) -> [f64; 2] {
    const RADIX: u32 = 2;
    const WIDTH: usize = 2;

    format!("{decimal:0WIDTH$b}")
            .chars()
            .flat_map(|c| c.to_digit(RADIX))
            .map(|n| n as f64)
            .collect::<Vec<f64>>()
            .try_into()
            .expect("Array wrong size compared to the iterator.")
}

#[allow(dead_code)]
pub fn base_counts_in_fasta(file: File) -> HashMap<char, u32> {
    let reader: BufReader<GzDecoder<File>> = BufReader::new(GzDecoder::new(file));
    let mut base_counts: HashMap<char, u32> = HashMap::new();

    for line in reader.lines().map(|l| l.unwrap()) {
        if (&line).starts_with(">") {
            continue;
        }

        for c in (&line).chars() {
            *base_counts.entry(c).or_insert(0) += 1;
        }
    }

    base_counts
}


#[derive(Debug, Hash, Copy, Clone, Eq, PartialEq)]
enum Nucleotide {
    A, 
    C, 
    G, 
    T, 
}


impl Distribution<Nucleotide> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Nucleotide {
        match rng.gen_range(0..=3) {
            0 => Nucleotide::A, 
            1 => Nucleotide::C, 
            2 => Nucleotide::G, 
            3 => Nucleotide::T, 
            _ => unreachable!(), 
        }
    }
}


fn str_to_nucleotides(s: &str) -> Vec<Nucleotide> {
    let choices: [Nucleotide; 4] = [Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T];
    let dist_n: WeightedIndex<u8> = WeightedIndex::new(&[1_u8, 1, 1, 1]).unwrap();
    let mut rng:ThreadRng = thread_rng();

    s.to_ascii_uppercase().chars().map(|base| match base {
        'A' => Nucleotide::A, 
        'C' => Nucleotide::C, 
        'G' => Nucleotide::G, 
        'T' => Nucleotide::T, 
        'N' => choices[dist_n.sample(&mut rng)],
        unk => panic!("Nucleotide base of unsupported type found in FASTA file: {}.", unk), 
    }).collect()
}

pub fn temp() -> () {
    println!("{:?}", str_to_nucleotides("nnnnatgccATTN"));
}


pub struct CGR {
    forward: Vec<f64>,
    backward: Vec<f64>,
}


pub fn compute_chaos_game_representation_from_fasta(file: File) -> CGR {
    let cgr_edges: HashMap<Nucleotide, [f64; 2]> = HashMap::from([
        (Nucleotide::A, [0.0, 1.0]),
    ]);
    let reader: BufReader<GzDecoder<File>> = BufReader::new(GzDecoder::new(file));

    let mut forward: Vec<f64> = Vec::new();
    let mut backward: Vec<f64> = Vec::new();

    for line in reader.lines().map(|l| l.unwrap()) {
        println!("{}", line);
    }
    
    CGR { forward, backward }
}


pub fn parse_fasta_file(filepath: &Path) -> std::io::Result<()> {
    let f: File = File::open(filepath).expect("Error in opening the FASTA file.");

    // let base_counts: HashMap<char, u32> = base_counts_in_fasta(f);
    // println!("{:?}", base_counts);
    // let n: u32 =  base_counts.iter().map(|(_, counts)| counts).sum();
    // println!("Sequence length = {}", n);

    compute_chaos_game_representation_from_fasta(f);

    Ok(())
}


