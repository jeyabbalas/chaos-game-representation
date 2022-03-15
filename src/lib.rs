use std::collections::HashMap;
use std::convert::TryInto;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use enum_iterator::IntoEnumIterator;
use rand::prelude::ThreadRng;
use rand::{
    distributions::{Distribution, Standard, WeightedIndex},
    Rng,
    thread_rng, 
};
use rev_lines::RevLines;


#[derive(Debug, Hash, Copy, Clone, IntoEnumIterator, PartialEq, Eq)]
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


pub struct ChaosGameRepresentation {
    forward: Vec<[f64;2]>,
    backward: Vec<[f64;2]>,
}


impl ChaosGameRepresentation {
    pub fn from_fasta_file(filename: &str) -> std::io::Result<ChaosGameRepresentation> {
        ChaosGameRepresentation::construct_chaos_game_representation(filename)
    }

    fn convert_decimal_to_binary(decimal: usize) -> [f64; 2] {
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

    fn compute_cgr_edges() -> HashMap<Nucleotide, [f64; 2]> {
        let mut cgr_edges: HashMap<Nucleotide, [f64; 2]> = HashMap::new();
    
        for (idx, base) in Nucleotide::into_enum_iter().enumerate() {
            cgr_edges.insert(base, ChaosGameRepresentation::convert_decimal_to_binary(idx));
        }
    
        cgr_edges
    }

    fn add(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
        let mut c: [f64; 2] = [0.0; 2];
        for (i, (a_val, b_val)) in a.iter().zip(&b).enumerate() {
            c[i] = a_val + b_val;
        }
        c
    }

    fn subtract(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
        let mut c: [f64; 2] = [0.0; 2];
        for (i, (a_val, b_val)) in a.iter().zip(&b).enumerate() {
            c[i] = a_val - b_val;
        }
        c
    }

    fn divide(a: [f64; 2], b: f64) -> [f64; 2] {
        let mut c: [f64; 2] = [0.0; 2];
        for (i, a_val) in a.iter().enumerate() {
            c[i] = a_val/b;
        }
        c
    }

    fn str_to_nucleotides(s: &str) -> Vec<Nucleotide> {
        let dna_bases: [Nucleotide; 4] = [Nucleotide::A, Nucleotide::C, Nucleotide::G, Nucleotide::T];
        let dist_n: WeightedIndex<u8> = WeightedIndex::new(&[1, 1, 1, 1]).unwrap();
        let mut rng: ThreadRng = thread_rng();
    
        s.to_ascii_uppercase().chars().map(|base| match base {
            'A' => Nucleotide::A, 
            'C' => Nucleotide::C, 
            'G' => Nucleotide::G, 
            'T' => Nucleotide::T, 
            'N' => dna_bases[dist_n.sample(&mut rng)],
            unk => panic!("Nucleotide base of unsupported type found in FASTA file: {}.", unk), 
        }).collect()
    }

    fn construct_chaos_game_representation(filename: &str) -> std::io::Result<ChaosGameRepresentation> {
        let cgr_edges: HashMap<Nucleotide, [f64; 2]> = ChaosGameRepresentation::compute_cgr_edges();
        let reader = BufReader::new(File::open(filename)?);
        let rev_lines = RevLines::new(BufReader::new(File::open(filename)?)).unwrap();

        let mut forward: Vec<[f64;2]> = Vec::new();
        let mut backward: Vec<[f64;2]> = Vec::new();

        let mut prev_point: [f64; 2] = [0.5; 2];

        let add = ChaosGameRepresentation::add;
        let subtract = ChaosGameRepresentation::subtract;
        let divide = ChaosGameRepresentation::divide;

        for line in reader.lines().map(|l| l.unwrap()) {
            if (&line).starts_with(">") {
                continue;
            }

            for base in ChaosGameRepresentation::str_to_nucleotides(&line) {
                prev_point = add(prev_point, 
                                 divide(subtract(cgr_edges[&base], prev_point), 2.0));
                forward.push(prev_point);
            }
        }

        for line in rev_lines {
            if (&line).len() == 0 || (&line).starts_with(">") {
                continue;
            }

            for base in ChaosGameRepresentation::str_to_nucleotides(&(&line).chars().rev().collect::<String>()) {
                prev_point = add(prev_point, 
                                 divide(subtract(cgr_edges[&base], prev_point), 2.0));
                backward.push(prev_point);
            }
        }
    
        Ok(ChaosGameRepresentation { forward, backward })
    }
}


pub fn run(filename_fasta: &str){
    let cgr = ChaosGameRepresentation::from_fasta_file(filename_fasta).expect("Error in parsing FASTA file.");

    println!("{}", &cgr.forward.len());
    println!("{}", &cgr.backward.len());
}


