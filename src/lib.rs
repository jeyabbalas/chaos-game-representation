pub mod fasta;

use std::collections::HashMap;
use std::fs::File;
use std::path::Path;
use std::io::prelude::*;
use std::io::BufReader;
use std::ops::{Add, Div, Sub};

use plotters::prelude::*;
use rand::prelude::ThreadRng;
use rand::{
    distributions::{
        Distribution, 
        Standard, 
        WeightedIndex
    },
    Rng,
    thread_rng, 
};
use rev_lines::RevLines;


#[derive(Debug, Hash, Copy, Clone, PartialEq, Eq)]
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
    name: String, 
    forward: Vec<Point<f64>>,
    backward: Vec<Point<f64>>,
}


impl ChaosGameRepresentation {
    pub fn from_fasta_file(filepath: &Path) -> Self {
        Self::construct_chaos_game_representation(filepath)
    }

    fn get_cgr_edges() -> HashMap<Nucleotide, Point<f64>> {
        use Nucleotide::*;

        let mut cgr_edges: HashMap<Nucleotide, Point<f64>> = HashMap::new();
    
        cgr_edges.insert(A, Point { x: 0.0, y: 0.0 } );
        cgr_edges.insert(C, Point { x: 0.0, y: 1.0 } );
        cgr_edges.insert(G, Point { x: 1.0, y: 0.0 } );
        cgr_edges.insert(T, Point { x: 1.0, y: 1.0 } );
    
        cgr_edges
    }

    fn str_to_nucleotides(s: &str) -> Vec<Nucleotide> {
        use Nucleotide::*;

        let dna_bases: [Nucleotide; 4] = [A, C, G, T];
        let dist_n: WeightedIndex<u8> = WeightedIndex::new(&[1, 1, 1, 1]).unwrap();
        let mut rng: ThreadRng = thread_rng();
    
        s.to_ascii_uppercase().chars().map(|base| match base {
            'A' => A, 
            'C' => C, 
            'G' => G, 
            'T' => T, 
            'N' => dna_bases[dist_n.sample(&mut rng)], 
            unk => panic!("Nucleotide base of unsupported type found in FASTA file: {}.", unk), 
        }).collect()
    }

    fn construct_chaos_game_representation(filepath: &Path) -> Self {
        let name = filepath.file_stem().unwrap().to_str().unwrap().to_string();

        let cgr_edges: HashMap<Nucleotide, Point<f64>> = ChaosGameRepresentation::get_cgr_edges();
        let reader = BufReader::new(File::open(filepath).expect("Error reading FASTA file."));
        let rev_lines = RevLines::new(BufReader::new(File::open(filepath).expect("Error reading FASTA file backwards."))).unwrap();

        let mut forward = Vec::new();
        let mut backward = Vec::new();

        let str_to_nucleotides = Self::str_to_nucleotides;

        let mut prev_point = Point { x: 0.5, y: 0.5 };

        for line in reader.lines().map(|l| l.unwrap()) {
            if (&line).starts_with(">") {
                continue;
            }

            for base in str_to_nucleotides(&line) {
                prev_point = prev_point + ((cgr_edges[&base] - prev_point) / 2.0);
                forward.push(prev_point);
            }
        }

        prev_point = Point { x: 0.5, y: 0.5 };

        for line in rev_lines {
            if (&line).len() == 0 || (&line).starts_with(">") {
                continue;
            }

            for base in str_to_nucleotides(&(&line).chars().rev().collect::<String>()) {
                prev_point = prev_point + ((cgr_edges[&base] - prev_point) / 2.0);
                backward.push(prev_point);
            }
        }
    
        Self { name, forward, backward }
    }
}


impl ChaosGameRepresentation {
    pub fn plot(&self, outdir: &Path, dimension: u32, margin: u32) -> Result<(), Box<dyn std::error::Error>> {
        let filepath_forward_cgr = outdir.join(format!("{}_forward.png", self.name));
        let filepath_backward_cgr = outdir.join(format!("{}_backward.png", self.name));

        self.plot_cgr(&filepath_forward_cgr,dimension, margin, true)?;
        self.plot_cgr(&filepath_backward_cgr, dimension, margin, false)?;

        Ok(())
    }

    fn plot_cgr(&self, filepath: &Path, dimension: u32, margin: u32, forward: bool) -> Result<(), Box<dyn std::error::Error>> {
        let cgr = if forward {
            &self.forward
        } else {
            &self.backward
        };

        let root = BitMapBackend::new(&filepath, (dimension+margin, dimension+margin)).into_drawing_area();
        root.fill(&WHITE)?;
        let root = root.margin(margin, margin, margin, margin);

        let mut chart = ChartBuilder::on(&root)
        .build_cartesian_2d(0.0..1.0, 1.0..0.0)?;

        chart.draw_series(cgr.iter().cloned().map(|point| Pixel::new((point.x, point.y), &BLACK)))?;
        Ok(())
    }
}


pub struct BufferedChaosGameRepresentation {
    fasta: fasta::Fasta,
}


impl BufferedChaosGameRepresentation {
    pub fn new(filepath: &Path) -> Self {
        let fasta = fasta::Fasta::new(filepath)
                .expect("Error landmarking FASTA file.");

        Self {
            fasta, 
        }
    }

    pub fn write_cgrs_to_hdf5() {
        // TODO
        // inputs: hdf5_filename, subset_seq_names &vec[1,2,3], chunks_shape
        // each sequence can be parallelized
    }

    pub fn plot_cgrs() {
        // TODO
    }
}

