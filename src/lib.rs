use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::Path;
use flate2::read::GzDecoder;
use std::str::from_utf8;


#[allow(dead_code)]
pub fn convert_decimal_to_binary(decimal: usize, width: usize) -> Vec<u32> {
    const RADIX: u32 = 2;

    format!("{decimal:0width$b}")
            .chars()
            .flat_map(|c| c.to_digit(RADIX))
            .collect()
}


pub fn decode_in_utf8(encoded_message: &[u8]) -> &str {
    match from_utf8(encoded_message) {
        Ok(decoded_message) => decoded_message,
        Err(e) => panic!("Invalid UTF-8 sequence found in FASTA file: {}.", e),
    }
}


pub fn parse_fasta_file() -> std::io::Result<()> {
    let fasta_filepath: &Path = Path::new("./data/egfr/egfr.fa.gz");
    let f: File = File::open(fasta_filepath).expect("Error reading FASTA file.");
    const CAPACITY: usize = 10*1000;
    let mut reader: BufReader<GzDecoder<File>> = BufReader::with_capacity(CAPACITY, GzDecoder::new(f));

    let mut start: bool = false;

    loop {
        let buffer: &[u8] = reader.fill_buf()?;
        let length: usize = buffer.len();

        if length == 0 {
            break;
        }

        let message = if !start {
            match buffer.iter().position(|&x| x == b'\n') {
                Some(idx) => {
                    start = true;
                    decode_in_utf8(&buffer[idx+1..])
                },
                None => {
                    continue;
                },
            }
        } else {
            decode_in_utf8(buffer)
        };

        println!("{}", message);

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


