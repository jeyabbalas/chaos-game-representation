use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::ops::Add;
use std::ops::Div;
use std::ops::Sub;
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


#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Point<T> {
    pub x: T, 
    pub y: T, 
}


impl<T: Add<Output = T>> Add for Point<T> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x + rhs.x, 
            y: self.y + rhs.y, 
        }
    }
}


impl<T: Sub<Output = T>> Sub for Point<T> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x - rhs.x, 
            y: self.y - rhs.y, 
        }
    }
}


impl<T: Div<Output = T> + Copy> Div<T> for Point<T> {
    type Output = Self;

    fn div(self, rhs: T) -> Self::Output {
        Self {
            x: self.x / rhs, 
            y: self.y / rhs, 
        }
    }
}


pub struct ChaosGameRepresentation {
    forward: Vec<Point<f64>>,
    backward: Vec<Point<f64>>,
}


impl ChaosGameRepresentation {
    pub fn from_fasta_file(filename: &str) -> ChaosGameRepresentation {
        ChaosGameRepresentation::construct_chaos_game_representation(filename)
    }

    fn compute_cgr_edges() -> HashMap<Nucleotide, Point<f64>> {
        let mut cgr_edges: HashMap<Nucleotide, Point<f64>> = HashMap::new();
    
        cgr_edges.insert(Nucleotide::A, Point { x: 0.0, y: 0.0} );
        cgr_edges.insert(Nucleotide::C, Point { x: 0.0, y: 1.0} );
        cgr_edges.insert(Nucleotide::G, Point { x: 1.0, y: 0.0} );
        cgr_edges.insert(Nucleotide::T, Point { x: 1.0, y: 1.0} );
    
        cgr_edges
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

    fn construct_chaos_game_representation(filename: &str) -> ChaosGameRepresentation {
        let cgr_edges: HashMap<Nucleotide, Point<f64>> = ChaosGameRepresentation::compute_cgr_edges();
        let reader = BufReader::new(File::open(filename).expect("Error reading FASTA file."));
        let rev_lines = RevLines::new(BufReader::new(File::open(filename).expect("Error reading FASTA file."))).unwrap();

        let mut forward = Vec::new();
        let mut backward = Vec::new();

        let mut prev_point = Point { x: 0.5, y: 0.5 };

        let str_to_nucleotides = ChaosGameRepresentation::str_to_nucleotides;

        for line in reader.lines().map(|l| l.unwrap()) {
            if (&line).starts_with(">") {
                continue;
            }

            for base in str_to_nucleotides(&line) {
                prev_point = prev_point + ((cgr_edges[&base] - prev_point) / 2.0);
                forward.push(prev_point);
            }
        }

        for line in rev_lines {
            if (&line).len() == 0 || (&line).starts_with(">") {
                continue;
            }

            for base in str_to_nucleotides(&(&line).chars().rev().collect::<String>()) {
                prev_point = prev_point + ((cgr_edges[&base] - prev_point) / 2.0);
                backward.push(prev_point);
            }
        }
    
        ChaosGameRepresentation { forward, backward }
    }
}


pub fn run(filename_fasta: &str){
    let cgr = ChaosGameRepresentation::from_fasta_file(filename_fasta);

    println!("{}", &cgr.forward.len());
    println!("{}", &cgr.backward.len());
}
