use std::io::prelude::*;
use std::fs::File;
use std::io::BufReader;


#[derive(Debug)]
pub struct FastaSequence {
    sequence_id: String, 
    start_byte: u64,
    end_byte: u64,
}


pub struct ForwardReader<R> {
    reader: BufReader<R>, 
    cursor_position: u64, 
    buf_size: u64, 
}


pub struct ReverseReader<R> {
    reader: BufReader<R>, 
    cursor_position: u64, 
    buf_size: u64, 
}



pub struct FastaReader {
    sequences: Vec<FastaSequence>, 
    //forward_reader: ForwardReader<R>, 
    //reverse_reader: ReverseReader<R>, 
}


impl FastaReader {
    pub fn new(reader: BufReader<File>) -> Result<FastaReader, Box<dyn std::error::Error>> {
        let mut sequences: Vec<FastaSequence> = Vec::new();
        
        let mut sequence_id = String::new();
        let mut start_byte = 0_u64;
        let mut end_byte: u64;

        let mut record = false;
        let mut cursor_position = 0_u64;
        let mut blank_line_counters = 0_u64;
        let mut unknown_id = 0_u32;

        for line in reader.lines().map(|l| l.unwrap()) {
            cursor_position += (line.len() + 1) as u64; // +1 for line feed character

            if line.len() == 0 {
                blank_line_counters += 1;
            }
            
            if line.starts_with(">") {
                if record {
                    // remove current header line chars + header line feed + 
                    // previous sequence end line feed + blank lines.
                    end_byte = cursor_position - (line.len() + 2) as u64 - blank_line_counters; 
                    sequences.push(FastaSequence {
                        sequence_id, 
                        start_byte, 
                        end_byte, 
                    });
                    blank_line_counters = 0;
                }
                
                sequence_id = (&line)[1..]
                        .split_whitespace()
                        .next()
                        .unwrap_or({
                            unknown_id = unknown_id + 1; 
                            &format!("Unknown-{}", unknown_id)
                        }).to_owned();
                start_byte = cursor_position;
                
                record = true;
            }
        }

        if record {
            // remove last line feed character + blank lines
            end_byte = cursor_position - 1 - blank_line_counters; 
            sequences.push(FastaSequence {
                sequence_id, 
                start_byte, 
                end_byte, 
            });
        }

        Ok(FastaReader {
            sequences
        })
    }
}


impl FastaReader {
    pub fn get_sequences(&self) -> &Vec<FastaSequence> {
        &self.sequences
    }
}



#[cfg(test)]
mod tests {
    use crate::utils::FastaReader;
    
    use std::io::prelude::*;
    use std::io::{BufReader, SeekFrom};
    use std::fs::File;
    use std::path::Path;

    #[test]
    pub fn parse_single_fasta() {
        let filepath = Path::new("./tests/data/single_fasta.fa");
        let reader = BufReader::new(File::open(filepath).expect("Error opening FASTA file."));
        let fasta_reader = FastaReader::new(reader).expect("FastaReader error");
        
        let sequences = fasta_reader.get_sequences();
        assert_eq!(sequences.len(), 1);
        assert_eq!(sequences[0].sequence_id, "sequenceID-001");
        
        let mut f = File::open("./tests/data/single_fasta.fa").expect("Error opening FASTA file.");
        let mut buffer = [0; 5];

        f.seek(SeekFrom::Start(sequences[0].start_byte)).expect("Error seeking.");
        f.read(&mut buffer[..]).expect("Error reading FASTA file");
        let seq = std::str::from_utf8(&buffer).expect("Invalid UTF-8 sequence.");
        assert_eq!(seq, "AAGTA");

        f.seek(SeekFrom::Start(sequences[0].end_byte - 5)).expect("Error seeking.");
        f.read(&mut buffer[..]).expect("Error reading FASTA file");
        let seq = std::str::from_utf8(&buffer).expect("Invalid UTF-8 sequence.");
        assert_eq!(seq, "CCCTG");
    }

    #[test]
    pub fn parse_anonymous_single_fasta() {
        let filepath = Path::new("./tests/data/anonymous_single_fasta.fa");
        let reader = BufReader::new(File::open(filepath).expect("Error opening FASTA file."));
        let fasta_reader = FastaReader::new(reader).expect("FastaReader error");
        
        let sequences = fasta_reader.get_sequences();
        assert_eq!(sequences.len(), 1);
        assert_eq!(sequences[0].sequence_id, "Unknown-1");
        
        let mut f = File::open("./tests/data/anonymous_single_fasta.fa").expect("Error opening FASTA file.");
        let mut buffer = [0; 5];

        f.seek(SeekFrom::Start(sequences[0].start_byte)).expect("Error seeking.");
        f.read(&mut buffer[..]).expect("Error reading FASTA file");
        let seq = std::str::from_utf8(&buffer).expect("Invalid UTF-8 sequence.");
        assert_eq!(seq, "AAGTA");

        f.seek(SeekFrom::Start(sequences[0].end_byte - 5)).expect("Error seeking.");
        f.read(&mut buffer[..]).expect("Error reading FASTA file");
        let seq = std::str::from_utf8(&buffer).expect("Invalid UTF-8 sequence.");
        assert_eq!(seq, "CCCTG");
    }

    #[test]
    pub fn parse_multi_fasta() {
        let filepath = Path::new("./tests/data/multi_fasta.fa");
        let reader = BufReader::new(File::open(filepath).expect("Error opening FASTA file."));
        let fasta_reader = FastaReader::new(reader).expect("FastaReader error");
        
        let sequences = fasta_reader.get_sequences();
        assert_eq!(sequences.len(), 3);
        assert_eq!(sequences[0].sequence_id, "sequenceID-001");
        assert_eq!(sequences[1].sequence_id, "sequenceID-002");
        assert_eq!(sequences[2].sequence_id, "sequenceID-003");
        
        let mut f = File::open("./tests/data/multi_fasta.fa").expect("Error opening FASTA file.");
        let mut buffer = [0; 5];

        f.seek(SeekFrom::Start(sequences[0].start_byte)).expect("Error seeking.");
        f.read(&mut buffer[..]).expect("Error reading FASTA file");
        let seq = std::str::from_utf8(&buffer).expect("Invalid UTF-8 sequence.");
        assert_eq!(seq, "AAGTA");

        f.seek(SeekFrom::Start(sequences[0].end_byte - 5)).expect("Error seeking.");
        f.read(&mut buffer[..]).expect("Error reading FASTA file");
        let seq = std::str::from_utf8(&buffer).expect("Invalid UTF-8 sequence.");
        assert_eq!(seq, "CCCTG");

        f.seek(SeekFrom::Start(sequences[1].start_byte)).expect("Error seeking.");
        f.read(&mut buffer[..]).expect("Error reading FASTA file");
        let seq = std::str::from_utf8(&buffer).expect("Invalid UTF-8 sequence.");
        assert_eq!(seq, "CAGTA");

        f.seek(SeekFrom::Start(sequences[1].end_byte - 5)).expect("Error seeking.");
        f.read(&mut buffer[..]).expect("Error reading FASTA file");
        let seq = std::str::from_utf8(&buffer).expect("Invalid UTF-8 sequence.");
        assert_eq!(seq, "TCGTC");

        f.seek(SeekFrom::Start(sequences[2].start_byte)).expect("Error seeking.");
        f.read(&mut buffer[..]).expect("Error reading FASTA file");
        let seq = std::str::from_utf8(&buffer).expect("Invalid UTF-8 sequence.");
        assert_eq!(seq, "CTTCA");

        f.seek(SeekFrom::Start(sequences[2].end_byte - 5)).expect("Error seeking.");
        f.read(&mut buffer[..]).expect("Error reading FASTA file");
        let seq = std::str::from_utf8(&buffer).expect("Invalid UTF-8 sequence.");
        assert_eq!(seq, "ATTTG");
    }
    
    #[test]
    pub fn parse_multi_fasta_extra_linefeeds() {
        let filepath = Path::new("./tests/data/multi_fasta_extra_linefeeds.fa");
        let reader = BufReader::new(File::open(filepath).expect("Error opening FASTA file."));
        let fasta_reader = FastaReader::new(reader).expect("FastaReader error");
        
        let sequences = fasta_reader.get_sequences();
        assert_eq!(sequences.len(), 3);
        assert_eq!(sequences[0].sequence_id, "sequenceID-001");
        assert_eq!(sequences[1].sequence_id, "sequenceID-002");
        assert_eq!(sequences[2].sequence_id, "sequenceID-003");
        
        let mut f = File::open("./tests/data/multi_fasta_extra_linefeeds.fa").expect("Error opening FASTA file.");
        let mut buffer = [0; 5];

        f.seek(SeekFrom::Start(sequences[0].start_byte)).expect("Error seeking.");
        f.read(&mut buffer[..]).expect("Error reading FASTA file");
        let seq = std::str::from_utf8(&buffer).expect("Invalid UTF-8 sequence.");
        assert_eq!(seq, "AAGTA");

        f.seek(SeekFrom::Start(sequences[0].end_byte - 5)).expect("Error seeking.");
        f.read(&mut buffer[..]).expect("Error reading FASTA file");
        let seq = std::str::from_utf8(&buffer).expect("Invalid UTF-8 sequence.");
        assert_eq!(seq, "CCCTG");

        f.seek(SeekFrom::Start(sequences[1].start_byte)).expect("Error seeking.");
        f.read(&mut buffer[..]).expect("Error reading FASTA file");
        let seq = std::str::from_utf8(&buffer).expect("Invalid UTF-8 sequence.");
        assert_eq!(seq, "CAGTA");

        f.seek(SeekFrom::Start(sequences[1].end_byte - 5)).expect("Error seeking.");
        f.read(&mut buffer[..]).expect("Error reading FASTA file");
        let seq = std::str::from_utf8(&buffer).expect("Invalid UTF-8 sequence.");
        assert_eq!(seq, "TCGTC");

        f.seek(SeekFrom::Start(sequences[2].start_byte)).expect("Error seeking.");
        f.read(&mut buffer[..]).expect("Error reading FASTA file");
        let seq = std::str::from_utf8(&buffer).expect("Invalid UTF-8 sequence.");
        assert_eq!(seq, "CTTCA");

        f.seek(SeekFrom::Start(sequences[2].end_byte - 5)).expect("Error seeking.");
        f.read(&mut buffer[..]).expect("Error reading FASTA file");
        let seq = std::str::from_utf8(&buffer).expect("Invalid UTF-8 sequence.");
        assert_eq!(seq, "ATTTG");
    }

    #[test]
    pub fn parse_anonymous_multi_fasta() {
        let filepath = Path::new("./tests/data/anonymous_multi_fasta.fa");
        let reader = BufReader::new(File::open(filepath).expect("Error opening FASTA file."));
        let fasta_reader = FastaReader::new(reader).expect("FastaReader error");
        
        let sequences = fasta_reader.get_sequences();
        assert_eq!(sequences.len(), 3);
        assert_eq!(sequences[0].sequence_id, "Unknown-1");
        assert_eq!(sequences[1].sequence_id, "Unknown-2");
        assert_eq!(sequences[2].sequence_id, "Unknown-3");
        
        let mut f = File::open("./tests/data/anonymous_multi_fasta.fa").expect("Error opening FASTA file.");
        let mut buffer = [0; 5];

        f.seek(SeekFrom::Start(sequences[0].start_byte)).expect("Error seeking.");
        f.read(&mut buffer[..]).expect("Error reading FASTA file");
        let seq = std::str::from_utf8(&buffer).expect("Invalid UTF-8 sequence.");
        assert_eq!(seq, "AAGTA");

        f.seek(SeekFrom::Start(sequences[0].end_byte - 5)).expect("Error seeking.");
        f.read(&mut buffer[..]).expect("Error reading FASTA file");
        let seq = std::str::from_utf8(&buffer).expect("Invalid UTF-8 sequence.");
        assert_eq!(seq, "CCCTG");

        f.seek(SeekFrom::Start(sequences[1].start_byte)).expect("Error seeking.");
        f.read(&mut buffer[..]).expect("Error reading FASTA file");
        let seq = std::str::from_utf8(&buffer).expect("Invalid UTF-8 sequence.");
        assert_eq!(seq, "CAGTA");

        f.seek(SeekFrom::Start(sequences[1].end_byte - 5)).expect("Error seeking.");
        f.read(&mut buffer[..]).expect("Error reading FASTA file");
        let seq = std::str::from_utf8(&buffer).expect("Invalid UTF-8 sequence.");
        assert_eq!(seq, "TCGTC");

        f.seek(SeekFrom::Start(sequences[2].start_byte)).expect("Error seeking.");
        f.read(&mut buffer[..]).expect("Error reading FASTA file");
        let seq = std::str::from_utf8(&buffer).expect("Invalid UTF-8 sequence.");
        assert_eq!(seq, "CTTCA");

        f.seek(SeekFrom::Start(sequences[2].end_byte - 5)).expect("Error seeking.");
        f.read(&mut buffer[..]).expect("Error reading FASTA file");
        let seq = std::str::from_utf8(&buffer).expect("Invalid UTF-8 sequence.");
        assert_eq!(seq, "ATTTG");
    }
}