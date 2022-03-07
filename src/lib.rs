use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::Path;
use flate2::read::GzDecoder;
use std::str::from_utf8;


#[allow(dead_code)]
pub fn convert_decimal_to_binary(decimal: usize, width: usize) -> Vec<u32> {
    const RADIX: u32 = 2;
    let binary: Vec<u32> = format!("{decimal:0width$b}")
                            .chars()
                            .flat_map(|c| c.to_digit(RADIX))
                            .collect();
    
    binary
}


pub fn parse_file() -> std::io::Result<()> {
    let fasta_filepath: &Path = Path::new("./data/egfr/egfr.fa.gz");
    let f: File = File::open(fasta_filepath).expect("Error reading FASTA file.");
    const CAPACITY: usize = 10*1000;
    let mut reader: BufReader<GzDecoder<File>> = BufReader::with_capacity(CAPACITY, GzDecoder::new(f));

    loop {
        let buffer: &[u8] = reader.fill_buf()?;
        let length: usize = buffer.len();

        if length == 0 {
            break;
        }

        let message: &str = match from_utf8(buffer) {
            Ok(buffer_decoded) => buffer_decoded,
            Err(e) => panic!("Invalid UTF-8 sequence found in FASTA file: {}.", e),
        };

        println!("{}", message.chars().count());

        reader.consume(length);
    }

    Ok(())
}


pub fn greet() -> () {
    let seq = "ATGCAT";

    for (i, c) in seq.chars().enumerate() {
        println!("Base# {i} is {c}", );
    }
}


