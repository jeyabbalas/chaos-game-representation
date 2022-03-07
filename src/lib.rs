use std::collections::HashSet;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::Path;
use flate2::read::GzDecoder;


#[allow(dead_code)]
pub fn convert_decimal_to_binary(decimal: usize, width: usize) -> Vec<u32> {
    const RADIX: u32 = 2;

    format!("{decimal:0width$b}")
            .chars()
            .flat_map(|c| c.to_digit(RADIX))
            .collect()
}


pub fn parse_fasta_file() -> std::io::Result<()> {
    let fasta_filepath: &Path = Path::new("./data/egfr/egfr.fa.gz");
    //let fasta_filepath: &Path = Path::new("./data/grch37/hg19.fa.gz");
    let f: File = File::open(fasta_filepath).expect("Error reading FASTA file.");
    let reader: BufReader<GzDecoder<File>> = BufReader::new(GzDecoder::new(f));

    let mut count: usize = 0;
    let mut unique_bases: HashSet<char> = HashSet::new();

    for line in reader.lines().map(|l| l.unwrap()) {
        if (&line).starts_with(">") {
            println!("{}", &line);
            continue;
        }

        unique_bases.extend((&line).chars());
        count += (&line).len();
    }

    println!("Number of bases = {}", count);
    println!("Unique bases = {:?}", unique_bases);

    Ok(())
}


pub fn greet() -> () {
    let seq = "ATGCAT";

    for (i, c) in seq.chars().enumerate() {
        println!("Base# {i} is {c}", );
    }
}


