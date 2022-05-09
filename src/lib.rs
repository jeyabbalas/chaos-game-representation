pub mod fasta;

use std::{
    collections::HashMap,
    env, 
    error::Error, 
    fs,
    io::{prelude::*, SeekFrom},
    ops::{Add, Div, Sub},
    path::Path
};

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

use crate::fasta::{Fasta, FastaSequence};


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


pub struct BufferedChaosGameRepresentation { 
    name: String, 
    fasta: Fasta,
}


impl BufferedChaosGameRepresentation {
    pub fn new(filepath: &Path) -> Self {
        let name = filepath.file_stem()
        .unwrap()
        .to_str()
        .unwrap()
        .to_string();
        
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
    
        cgr_edges.insert(A, Point { x: 0.0, y: 1.0 } );
        cgr_edges.insert(C, Point { x: 1.0, y: 1.0 } );
        cgr_edges.insert(G, Point { x: 0.0, y: 0.0 } );
        cgr_edges.insert(T, Point { x: 1.0, y: 0.0 } );
    
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

    pub fn build_cgrs_and_write_to_binary_file(&self, output_dir: &Path) -> std::io::Result<()> {
            let sequence_ids = self.get_fasta().get_sequence_ids();
            self.build_select_cgrs_and_write_to_binary_file(output_dir, sequence_ids)
    }
    
    pub fn build_select_cgrs_and_write_to_binary_file(&self, 
        output_dir: &Path, 
        sequence_ids: Vec<&str>) -> std::io::Result<()> {
        fs::create_dir_all(output_dir)?;

        let cgr_edges: HashMap<Nucleotide, Point<f64>> = Self::get_cgr_edges();
        let str_to_nucleotides = Self::str_to_nucleotides;
        
        let mut sequences = Vec::<&FastaSequence>::with_capacity(sequence_ids.len());
        for sequence_id in sequence_ids.clone() {
            let fasta_seq = self.get_fasta()
                .get_sequence_by_id(sequence_id)
                .expect(&format!("Entered sequence ID: {sequence_id} not found"));
            sequences.push(fasta_seq);
        }
        
        for (sequence_id, sequence) in sequence_ids.into_iter().zip(sequences.into_iter()) {
            let sequence_len = sequence.get_sequence_length() as usize;
            println!("Writing sequence {} of length {} to binary file.", sequence_id, sequence_len);

            let sub_directory = output_dir.join(sequence_id);
            fs::create_dir_all(&sub_directory)?;

            // binary output files
            let forward_cgr_filename = sub_directory.join("forward_cgr.bin");
            let backward_cgr_filename = sub_directory.join("backward_cgr.bin");
            let forward_backward_cgr_filename = sub_directory.join("forward_backward_cgr.bin");
            
            let forward_reader = sequence.build_forward_reader()
                .expect("Error building forward reader");
            let backward_reader = sequence.build_reverse_reader()
                .expect("Error building reverse reader");
            
            let chunk_length = if sequence_len > 1_000_000 { 1_000_000 } else { sequence_len };
            const BYTES_IN_F64: usize = 8;
            const CGR_DIMENSION: usize = 2;
            {
                // write sequence and calculate and write forward cgr
                let mut sequence_bin = fs::File::create(sub_directory.join("sequence.bin"))?;
                let mut forward_cgr_bin = fs::File::create(&forward_cgr_filename)?;

                let mut prev_point = Point { x: 0.5, y: 0.5 };
                let mut sequence_parsed = 0;
                let mut sequence_bytes = Vec::<u8>::with_capacity(chunk_length);
                let mut sequence_buffer = Vec::<Nucleotide>::with_capacity(chunk_length);
                for line in forward_reader {
                    let mut line_bytes: Vec<u8> = (&line).as_bytes().to_vec();
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

                    let mut forward_cgr_bytes = Vec::<u8>::with_capacity(forward_cgr.len() * CGR_DIMENSION * BYTES_IN_F64);

                    for point in forward_cgr.iter() {
                        forward_cgr_bytes.extend_from_slice(&point.x.to_be_bytes());
                        forward_cgr_bytes.extend_from_slice(&point.y.to_be_bytes());
                    }

                    sequence_bin.write_all(&sequence_bytes[..])?;
                    forward_cgr_bin.write_all(&forward_cgr_bytes[..])?;

                    sequence_bytes = Vec::<u8>::with_capacity(chunk_length);
                    sequence_buffer = Vec::<Nucleotide>::with_capacity(chunk_length);
                }
            }

            {
                // calculate and write backward cgr
                let mut backward_cgr_bin = fs::File::create(&backward_cgr_filename)?;

                let mut prev_point = Point { x: 0.5, y: 0.5 };
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

                    let mut backward_cgr_bytes = Vec::<u8>::with_capacity(backward_cgr.len() * CGR_DIMENSION * BYTES_IN_F64);

                    for point in backward_cgr.iter() {
                        backward_cgr_bytes.extend_from_slice(&point.y.to_be_bytes());
                        backward_cgr_bytes.extend_from_slice(&point.x.to_be_bytes());
                    }

                    backward_cgr_bin.write_all(&backward_cgr_bytes[..])?;

                    sequence_buffer = Vec::<Nucleotide>::with_capacity(chunk_length);
                }
            }

            {
                let mut forward_cgr_bin = fs::File::open(&forward_cgr_filename)?;
                let mut backward_cgr_bin = fs::File::open(&backward_cgr_filename)?;
                let mut forward_backward_cgr_bin = fs::File::create(&forward_backward_cgr_filename)?;

                let mut sequence_parsed = 0_usize;

                while (sequence_len - sequence_parsed) > 0 {
                    let buffer_size = if (sequence_len - sequence_parsed) > chunk_length { chunk_length } else { sequence_len - sequence_parsed };

                    let mut forward_buffer = vec![0_u8; buffer_size * CGR_DIMENSION * BYTES_IN_F64];
                    let mut backward_buffer = vec![0_u8; buffer_size * CGR_DIMENSION * BYTES_IN_F64];

                    forward_cgr_bin.seek(SeekFrom::Start((sequence_parsed * CGR_DIMENSION * BYTES_IN_F64) as u64))?;
                    backward_cgr_bin.seek(SeekFrom::End(-(((sequence_parsed + buffer_size) * CGR_DIMENSION * BYTES_IN_F64) as i64)))?;

                    forward_cgr_bin.read_exact(&mut forward_buffer)?;
                    backward_cgr_bin.read_exact(&mut backward_buffer)?;

                    let mut forward_backward_cgr_bytes = Vec::<u8>::with_capacity(buffer_size * CGR_DIMENSION * BYTES_IN_F64 * 2);

                    for i in 0..buffer_size {
                        let forward_ptr = i * CGR_DIMENSION * BYTES_IN_F64;
                        let backward_ptr = (buffer_size - 1 - i) * CGR_DIMENSION * BYTES_IN_F64;

                        forward_backward_cgr_bytes.extend_from_slice(&forward_buffer[forward_ptr..(forward_ptr + (CGR_DIMENSION * BYTES_IN_F64))]);
                        forward_backward_cgr_bytes.extend_from_slice(&backward_buffer[backward_ptr..(backward_ptr + (CGR_DIMENSION * BYTES_IN_F64))]);
                    }

                    forward_backward_cgr_bin.write_all(&forward_backward_cgr_bytes[..])?;

                    sequence_parsed += buffer_size;
                }

                fs::remove_file(&forward_cgr_filename)?;
                fs::remove_file(&backward_cgr_filename)?;
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

    pub fn get_sequence_ids(&self) -> Vec<&str> {
        self.get_fasta().get_sequence_ids()
    }
}

pub struct Config {
    filename_fasta: String, 
    output_dir: String, 
    sequence_ids: Vec<String>, 
}

impl Config {
    pub fn new(mut args: env::Args) -> Result<Self, &'static str> {
        args.next();

        let filename_fasta = match args.next() {
            Some(f) => {f}, 
            None => {return Err("Error: input fasta file not specified.");}, 
        };

        let output_dir = match args.next() {
            Some(d) => {d}, 
            None => {return Err("Error: output directory not specified.");}, 
        };

        let sequence_ids: Vec<String> = match args.next() {
            Some(ids) => {
                ids.split(",").into_iter().map(|s| s.to_string()).collect()
            },
            None => {return Err("Error: list of fasta sequences to write to binary file not specified.");}, 
        };

        println!("Input file: {filename_fasta}");
        println!("Output directory: {output_dir}");
        println!("Sequences to write to binary file: {}.\n", sequence_ids.join(", "));

        Ok(Self{
            filename_fasta, 
            output_dir, 
            sequence_ids, 
        })
    }
}

pub fn run(config: Config) -> Result<(), Box<dyn Error>> {
    let filepath = Path::new(&config.filename_fasta);
    let output_dir = Path::new(&config.output_dir);

    let cgr_buf = BufferedChaosGameRepresentation::new(filepath);
    let sequence_ids: Vec<&str> = config.sequence_ids.iter().map(|s| s.as_str()).collect();
    cgr_buf.build_select_cgrs_and_write_to_binary_file(output_dir, sequence_ids)?;

    Ok(())
}
