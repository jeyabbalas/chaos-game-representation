pub mod fasta;

use std::collections::HashMap;
use std::path::Path;
use std::ops::{Add, Div, Sub};

use hdf5::{self, types::FixedAscii};
use ndarray::{Array1, arr2};
use plotters::prelude::*;
use rand::{
    prelude::ThreadRng, 
    distributions::{
        Distribution, 
        Standard, 
        WeightedIndex
    },
    Rng,
    thread_rng, 
};

use crate::fasta::Fasta;


#[derive(Debug, Hash, Copy, Clone, PartialEq, Eq)]
enum Nucleotide {
    A, 
    C, 
    G, 
    T, 
}


impl Distribution<Nucleotide> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Nucleotide {
        use Nucleotide::*;

        match rng.gen_range(0..=3) {
            0 => A, 
            1 => C, 
            2 => G, 
            3 => T, 
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
        let mut forward = Vec::<Point<f64>>::new();
        let mut backward = Vec::<Point<f64>>::new();

        let cgr_edges: HashMap<Nucleotide, Point<f64>> = Self::get_cgr_edges();

        let fasta = Fasta::new(filepath)
                .expect("Error landmarking FASTA file.");
        let sequences = fasta.get_sequences();
        let forward_reader = sequences[0].build_forward_reader()
                .expect("Error building forward reader");
        let backward_reader = sequences[0].build_reverse_reader()
                .expect("Error building reverse reader");

        let str_to_nucleotides = Self::str_to_nucleotides;

        let mut prev_point = Point { x: 0.5, y: 0.5 };

        for line in forward_reader {
            if (&line).starts_with(">") {
                continue;
            }

            for base in str_to_nucleotides(&line) {
                prev_point = prev_point + ((cgr_edges[&base] - prev_point) / 2.0);
                forward.push(prev_point);
            }
        }

        prev_point = Point { x: 0.5, y: 0.5 };

        for line in backward_reader {
            if (&line).starts_with(">") {
                continue;
            }

            for base in str_to_nucleotides(&line) {
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
    name: String, 
    fasta: Fasta,
}


impl BufferedChaosGameRepresentation {
    pub fn new(filepath: &Path) -> Self {
        let name = filepath.file_stem().unwrap().to_str().unwrap().to_string();
        let fasta = Fasta::new(filepath)
                .expect("Error landmarking FASTA file.");

        Self {
            name, 
            fasta, 
        }
    }
}


impl BufferedChaosGameRepresentation {
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
    
    pub fn build_cgrs_and_write_to_hdf5(&self, filename: &Path, chunk_length: usize) -> hdf5::Result<(), hdf5::Error> {
        let file = hdf5::File::create(filename)?;

        let cgr_edges: HashMap<Nucleotide, Point<f64>> = Self::get_cgr_edges();
        let str_to_nucleotides = Self::str_to_nucleotides;

        for sequence in self.get_fasta().get_sequences() {
            let name = sequence.get_sequence_id();
            let group = file.create_group(name)?;

            let len = sequence.get_sequence_length() as usize;

            // sequence
            let ds_sequence = group.new_dataset::<FixedAscii<1>>()
                .chunk(chunk_length)
                .shape(len)
                .deflate(5)
                .create("sequence")?;

            // forward cgr
            let ds_forward_cgr = group.new_dataset::<f64>()
                .chunk((chunk_length, 2))
                .shape((len, 2))
                .deflate(5)
                .create("forward_cgr")?;
            
            // backward cgr
            let ds_backward_cgr = group.new_dataset::<f64>()
                .chunk((chunk_length, 2))
                .shape((len, 2))
                .deflate(5)
                .create("backward_cgr")?;
            
            let forward_reader = sequence.build_forward_reader()
                .expect("Error building forward reader");
            let backward_reader = sequence.build_reverse_reader()
                .expect("Error building reverse reader");
            
            
            let mut prev_point = Point { x: 0.5, y: 0.5 };
            let mut idx = 0;
            for line in forward_reader {
                let mut forward_cgr = Vec::<Point<f64>>::new();
                
                if (&line).starts_with(">") {
                    continue;
                }
        
                for base in str_to_nucleotides(&line) {
                    prev_point = prev_point + ((cgr_edges[&base] - prev_point) / 2.0);
                    forward_cgr.push(prev_point);
                }

                let sequence_bytes: Vec<FixedAscii<1>> = (&line).to_ascii_uppercase()
                    .chars()
                    .map(|base| FixedAscii::<1>::from_ascii(&[base as u8]).unwrap())
                    .collect();
                let seq_arr = Array1::from(sequence_bytes);
                ds_sequence.write_slice(&seq_arr, idx..(idx + seq_arr.len()))?;

                let forward_coords: Vec<[f64; 2]> = forward_cgr.iter()
                    .map(|point| [point.x, point.y])
                    .collect();
                let forward_coords: &[[f64; 2]] = &forward_coords;
                let forward_arr = arr2(forward_coords);
                ds_forward_cgr.write_slice(&forward_arr, (idx..(idx + seq_arr.len()), ..) )?;

                idx += seq_arr.len();
            }

            let mut prev_point = Point { x: 0.5, y: 0.5 };
            let mut idx = 0;
            for line in backward_reader {
                let mut backward_cgr = Vec::<Point<f64>>::new();

                if (&line).starts_with(">") {
                    continue;
                }
    
                for base in str_to_nucleotides(&line) {
                    prev_point = prev_point + ((cgr_edges[&base] - prev_point) / 2.0);
                    backward_cgr.push(prev_point);
                }

                let backward_coords: Vec<[f64; 2]> = backward_cgr.iter()
                    .map(|point| [point.x, point.y])
                    .collect();
                let backward_coords: &[[f64; 2]] = &backward_coords;
                let backward_arr = arr2(backward_coords);
                let seq_length = backward_arr.shape()[0];
                ds_backward_cgr.write_slice(&backward_arr, (idx..(idx + seq_length), ..))?;

                idx += seq_length;
            }
        }

        Ok(())
    }
}


impl BufferedChaosGameRepresentation {
    pub fn get_name(&self) -> String {
        self.name.to_string()
    }

    pub fn get_fasta(&self) -> &Fasta {
        &self.fasta
    }

    pub fn get_sequence_ids(&self) -> Vec<&String> {
        self.get_fasta().get_sequence_ids()
    }
}

