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
                .expect("Error while indexing FASTA file.");

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

    pub fn build_cgrs_and_write_to_hdf5(&self, 
        filename: &Path, 
        chunk_length: usize) -> hdf5::Result<(), hdf5::Error> {
            let select_ids = self.get_fasta().get_sequence_ids();
            self.build_select_cgrs_and_write_to_hdf5(filename, chunk_length, select_ids)
    }
    
    pub fn build_select_cgrs_and_write_to_hdf5(&self, 
        filename: &Path, 
        chunk_length: usize, 
        select_ids: Vec<&str>) -> hdf5::Result<(), hdf5::Error> {
        use Nucleotide::*;
        let file = hdf5::File::create(filename)?;

        let cgr_edges: HashMap<Nucleotide, Point<f64>> = Self::get_cgr_edges();
        let str_to_nucleotides = Self::str_to_nucleotides;
        
        for sequence_id in select_ids {
            let sequence = self.get_fasta().get_sequence_by_id(sequence_id)
                .expect("Sequence ID not found.");
            let sequence_len = sequence.get_sequence_length() as usize;
            println!("Writing sequence {} of length {} to HDF5 file.", sequence_id, sequence_len);

            let group = file.create_group(sequence_id)?;

            // sequence
            let ds_sequence = group.new_dataset::<FixedAscii<1>>()
                .chunk(chunk_length)
                .shape(sequence_len)
                .deflate(5)
                .create("sequence")?;

            // forward cgr
            let ds_forward_cgr = group.new_dataset::<f64>()
                .chunk((chunk_length, 2))
                .shape((sequence_len, 2))
                .deflate(5)
                .create("forward_cgr")?;
            
            // backward cgr
            let ds_backward_cgr = group.new_dataset::<f64>()
                .chunk((chunk_length, 2))
                .shape((sequence_len, 2))
                .deflate(5)
                .create("backward_cgr")?;
            
            let forward_reader = sequence.build_forward_reader()
                .expect("Error building forward reader");
            let backward_reader = sequence.build_reverse_reader()
                .expect("Error building reverse reader");
            
            // write sequence and calculate and write forward cgr
            let mut prev_point = Point { x: 0.5, y: 0.5 };
            let mut idx = 0;
            let mut sequence_parsed = 0;
            let mut sequence_bytes = Vec::<FixedAscii<1>>::with_capacity(chunk_length);
            let mut sequence_buffer = Vec::<Nucleotide>::with_capacity(chunk_length);
            for line in forward_reader {
                let mut line_bytes: Vec<FixedAscii<1>> = (&line).to_ascii_uppercase()
                    .chars()
                    .map(|base| FixedAscii::<1>::from_ascii(&[base as u8]).unwrap())
                    .collect();
                sequence_bytes.append(&mut line_bytes);
                
                let mut line_nucleotides = str_to_nucleotides(&line);
                let line_len = line_nucleotides.len();
                sequence_buffer.append(&mut line_nucleotides);
                sequence_parsed += line_len;

                if sequence_buffer.len() < sequence_buffer.capacity() && 
                sequence_parsed < sequence_len {
                    continue;
                }

                let mut forward_cgr = Vec::<Point<f64>>::with_capacity(sequence_buffer.len());
        
                for base in sequence_buffer {
                    prev_point = prev_point + ((cgr_edges[&base] - prev_point) / 2.0);
                    forward_cgr.push(prev_point);
                }

                let seq_arr = Array1::from(sequence_bytes);
                ds_sequence.write_slice(&seq_arr, idx..(idx + seq_arr.len()))?;

                let forward_coords: Vec<[f64; 2]> = forward_cgr.iter()
                    .map(|point| [point.x, point.y])
                    .collect();
                let forward_coords: &[[f64; 2]] = &forward_coords;
                let forward_arr = arr2(forward_coords);
                ds_forward_cgr.write_slice(&forward_arr, (idx..(idx + seq_arr.len()), ..) )?;
                
                sequence_buffer = Vec::<Nucleotide>::with_capacity(chunk_length);
                sequence_bytes = Vec::<FixedAscii<1>>::with_capacity(chunk_length);
                idx += seq_arr.len();
            }

            let attr = ds_forward_cgr.new_attr::<f64>()
                .shape([2])
                .create("A")?;
            attr.write(&[cgr_edges[&A].x, cgr_edges[&A].y])?;
            
            let attr = ds_forward_cgr.new_attr::<f64>()
                .shape([2])
                .create("T")?;
            attr.write(&[cgr_edges[&T].x, cgr_edges[&T].y])?;
            
            let attr = ds_forward_cgr.new_attr::<f64>()
                .shape([2])
                .create("G")?;
            attr.write(&[cgr_edges[&G].x, cgr_edges[&G].y])?;
            
            let attr = ds_forward_cgr.new_attr::<f64>()
                .shape([2])
                .create("C")?;
            attr.write(&[cgr_edges[&C].x, cgr_edges[&C].y])?;

            // calculate and write backward cgr
            let mut prev_point = Point { x: 0.5, y: 0.5 };
            let mut idx = sequence_len;
            let mut sequence_parsed = 0;
            let mut sequence_buffer = Vec::<Nucleotide>::with_capacity(chunk_length);
            for line in backward_reader {
                let mut line_nucleotides = str_to_nucleotides(&line);
                let line_len = line_nucleotides.len();
                sequence_buffer.append(&mut line_nucleotides);
                sequence_parsed += line_len;

                if sequence_buffer.len() < sequence_buffer.capacity() && 
                sequence_parsed < sequence_len {
                    continue;
                }
                
                let mut backward_cgr = Vec::<Point<f64>>::new();
    
                for base in sequence_buffer {
                    prev_point = prev_point + ((cgr_edges[&base] - prev_point) / 2.0);
                    backward_cgr.push(prev_point);
                }

                let backward_coords: Vec<[f64; 2]> = backward_cgr.iter()
                    .map(|point| [point.x, point.y])
                    .rev()
                    .collect();
                let backward_coords: &[[f64; 2]] = &backward_coords;
                let backward_arr = arr2(backward_coords);
                let segment_length = backward_arr.shape()[0];
                ds_backward_cgr.write_slice(&backward_arr, ((idx - segment_length)..idx, ..))?;

                sequence_buffer = Vec::<Nucleotide>::with_capacity(chunk_length);
                idx -= segment_length;
            }

            let attr = ds_backward_cgr.new_attr::<f64>()
                .shape([2])
                .create("A")?;
            attr.write(&[cgr_edges[&A].x, cgr_edges[&A].y])?;
            
            let attr = ds_backward_cgr.new_attr::<f64>()
                .shape([2])
                .create("T")?;
            attr.write(&[cgr_edges[&T].x, cgr_edges[&T].y])?;
            
            let attr = ds_backward_cgr.new_attr::<f64>()
                .shape([2])
                .create("G")?;
            attr.write(&[cgr_edges[&G].x, cgr_edges[&G].y])?;
            
            let attr = ds_backward_cgr.new_attr::<f64>()
                .shape([2])
                .create("C")?;
            attr.write(&[cgr_edges[&C].x, cgr_edges[&C].y])?;
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

    pub fn get_sequence_ids(&self) -> Vec<&str> {
        self.get_fasta().get_sequence_ids()
    }
}

