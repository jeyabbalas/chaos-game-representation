use std::cmp::min;
use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader, Result, SeekFrom};
use std::path::{Path, PathBuf};


static DEFAULT_SIZE: usize = 4096;

static LF_BYTE: u8 = '\n' as u8;
static CR_BYTE: u8 = '\r' as u8;


pub struct ForwardFastaReader<R> {
    reader: BufReader<R>, 
    end_byte: u64,
    cursor_position: u64, 
    buffer_size: u64, 
}


impl<R:Seek+Read> ForwardFastaReader<R> {
    pub fn new(reader: BufReader<R>, fasta_sequence: &FastaSequence) -> Result<ForwardFastaReader<R>> {
        ForwardFastaReader::with_capacity(DEFAULT_SIZE, reader, 
            fasta_sequence.get_start_byte(), fasta_sequence.get_end_byte())
    }

    pub fn with_capacity(capacity: usize, mut reader: BufReader<R>, 
                         start_byte: u64, end_byte: u64) -> Result<ForwardFastaReader<R>> {
        if start_byte > end_byte {
            panic!("Start position of read cannot come after the end position of the sequence.");
        }

        let cursor_position = reader.seek(SeekFrom::Start(start_byte))?;

        Ok(ForwardFastaReader {
            reader, 
            end_byte, 
            cursor_position, 
            buffer_size: capacity as u64, 
        })
    }

    fn read_to_buffer(&mut self, size: u64) -> Result<Vec<u8>> {
        let mut buffer = vec![0; size as usize]; 
        self.reader.read_exact(&mut buffer[0..(size as usize)])?;
        self.cursor_position += size;

        Ok(buffer)
    }
}


impl<R:Read+Seek> Iterator for ForwardFastaReader<R> {
    type Item = String;
    
    fn next(&mut self) -> Option<String> {
        let mut result: Vec<u8> = Vec::new();

        'outer: loop {
            if self.cursor_position >= self.end_byte {
                if result.len() > 0 {
                    break;
                }
                return None;
            }

            let size = min(self.buffer_size, self.end_byte - self.cursor_position);

            match self.read_to_buffer(size) {
                Ok(buffer) => {
                    for (idx, c) in (&buffer).iter().enumerate() {
                        if *c == LF_BYTE || *c == CR_BYTE {
                            let mut offset = size - 1 - idx as u64;

                            // handling CRLF combos
                            if buffer[idx] == CR_BYTE && (idx as u64) < (size - 1) && buffer[idx] == LF_BYTE {
                                offset += 1;
                            }

                            match self.reader.seek(SeekFrom::Current(-(offset as i64))) {
                                Ok(_) => {
                                    self.cursor_position -= offset;
                                    break 'outer;
                                },
                                Err(_) => return None,
                            }
                        } else {
                            result.push(c.clone());
                        }
                    }
                },
                Err(_) => return None,
            }
        }

        Some(String::from_utf8(result).unwrap())
    }
}


pub struct ReverseFastaReader<R> {
    reader: BufReader<R>, 
    start_byte: u64,
    cursor_position: u64, 
    buffer_size: u64, 
}


impl<R:Seek+Read> ReverseFastaReader<R> {
    pub fn new(reader: BufReader<R>, fasta_sequence: &FastaSequence) -> Result<ReverseFastaReader<R>> {
        ReverseFastaReader::with_capacity(DEFAULT_SIZE, reader, 
            fasta_sequence.get_start_byte(), fasta_sequence.get_end_byte())
    }

    pub fn with_capacity(capacity: usize, mut reader: BufReader<R>, 
                         start_byte: u64, end_byte: u64) -> Result<ReverseFastaReader<R>> {
        if start_byte > end_byte {
            panic!("Start position of read cannot come after the end position of the sequence.");
        }

        let cursor_position = reader.seek(SeekFrom::Start(end_byte))?;

        Ok(ReverseFastaReader {
            reader, 
            start_byte, 
            cursor_position,  
            buffer_size: capacity as u64, 
        })
    }

    fn read_to_buffer(&mut self, size: u64) -> Result<Vec<u8>> {
        let mut buffer = vec![0; size as usize]; 
        let offset = -(size as i64);

        self.reader.seek(SeekFrom::Current(offset))?;
        self.reader.read_exact(&mut buffer[0..(size as usize)])?;
        self.reader.seek(SeekFrom::Current(offset))?;

        self.cursor_position -= size;

        Ok(buffer)
    }
}


impl<R:Read+Seek> Iterator for ReverseFastaReader<R> {
    type Item = String;
    
    fn next(&mut self) -> Option<String> {
        let mut result: Vec<u8> = Vec::new();

        'outer: loop {
            if self.cursor_position <= self.start_byte {
                if result.len() > 0 {
                    break;
                }
                return None;
            }

            let size = min(self.buffer_size, self.cursor_position - self.start_byte);

            match self.read_to_buffer(size) {
                Ok(buffer) => {
                    for (idx, c) in (&buffer).iter().enumerate().rev() {
                        if *c == LF_BYTE {
                            let mut offset = idx as u64;

                            // handling CRLF combos
                            if idx > 1 && buffer[idx-1] == CR_BYTE {
                                offset -= 1;
                            }

                            match self.reader.seek(SeekFrom::Current(offset as i64)) {
                                Ok(_) => {
                                    self.cursor_position += offset;
                                    break 'outer;
                                },
                                Err(_) => return None,
                            }
                        } else {
                            result.push(c.clone());
                        }
                    }
                },
                Err(_) => return None,
            }
        }

        Some(String::from_utf8(result).unwrap())
    }
}


pub struct FastaSequence {
    filepath: PathBuf, 
    sequence_id: String, 
    start_byte: u64, 
    end_byte: u64, 
    sequence_length: u64, 
}


impl FastaSequence {
    pub fn get_sequence_id(&self) -> &String {
        &self.sequence_id
    }

    pub fn get_start_byte(&self) -> u64 {
        self.start_byte
    }

    pub fn get_end_byte(&self) -> u64 {
        self.end_byte
    }

    pub fn get_sequence_length(&self) -> u64 {
        self.sequence_length
    }
}


impl FastaSequence {
    pub fn build_forward_reader(&self) -> Result<ForwardFastaReader<File>> {
        let reader = BufReader::new(
            File::open(&self.filepath)
            .expect("Error opening FASTA file when building forward reader."));

        ForwardFastaReader::new(reader, self)
    }

    pub fn build_reverse_reader(&self) -> Result<ReverseFastaReader<File>> {
        let reader = BufReader::new(
            File::open(&self.filepath)
            .expect("Error opening FASTA file when building reverse reader."));

        ReverseFastaReader::new(reader, self)
    }
}


pub struct Fasta {
    sequences: Vec<FastaSequence>, 
}


impl Fasta {
    pub fn new(fasta_filepath: &Path) -> Result<Fasta> {
        let reader = BufReader::new(
            File::open(fasta_filepath)
            .expect("Error opening FASTA file."));
        let mut sequences: Vec<FastaSequence> = Vec::new();
        
        let mut sequence_id = String::new();
        let mut start_byte = 0_u64;
        let mut end_byte: u64;
        let mut sequence_length = 0_u64;

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
                    // remove— current header line chars + header line feed + 
                    // previous sequence end line feed + blank lines.
                    end_byte = cursor_position - (line.len() + 2) as u64 - blank_line_counters; 
                    
                    sequences.push(FastaSequence {
                        filepath: fasta_filepath.to_path_buf(), 
                        sequence_id, 
                        start_byte, 
                        end_byte, 
                        sequence_length, 
                    });

                    sequence_length = 0;
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
            } else {
                sequence_length += line.len() as u64;
            }
        }

        if record {
            // remove last line feed character + blank lines
            end_byte = cursor_position - 1 - blank_line_counters; 

            sequences.push(FastaSequence {
                filepath: fasta_filepath.to_path_buf(), 
                sequence_id, 
                start_byte, 
                end_byte, 
                sequence_length, 
            });
        }

        Ok(Fasta {
            sequences, 
        })
    }
}


impl Fasta {
    pub fn get_sequences(&self) -> &Vec<FastaSequence> {
        return &self.sequences
    }

    pub fn get_sequence_ids(&self) -> Vec<&String> {
        self.get_sequences().iter().map(|seq| seq.get_sequence_id()).collect()
    }
}


#[cfg(test)]
mod tests {
    use crate::fasta::Fasta;
    use std::path::Path;

    #[test]
    pub fn parse_single_fasta() {
        let filepath = Path::new("./tests/data/single_fasta.fa");
        let fasta_reader = Fasta::new(filepath).expect("Error landmarking FASTA file.");
        
        let sequences = fasta_reader.get_sequences();
        assert_eq!(sequences.len(), 1);
        assert_eq!(sequences[0].get_sequence_id(), "sequenceID-001");
        assert_eq!(sequences[0].get_sequence_length(), 106);
        
        let mut forward_reader = sequences[0].build_forward_reader().expect("Error building forward reader");

        assert_eq!(forward_reader.next(), Some("AAGTAGGAATAATATCTTATCATTATAGATAAAAACCTTCTGAATTTGCTTAGTGTGTAT".to_string()));
        assert_eq!(forward_reader.next(), Some("ACGACTAGACATATATCAGCTCGCCGATTATTTGGATTATTCCCTG".to_string()));
        assert_eq!(forward_reader.next(), None);

        let mut reverse_reader = sequences[0].build_reverse_reader().expect("Error building reverse reader");

        assert_eq!(reverse_reader.next(), Some("GTCCCTTATTAGGTTTATTAGCCGCTCGACTATATACAGATCAGCA".to_string()));
        assert_eq!(reverse_reader.next(), Some("TATGTGTGATTCGTTTAAGTCTTCCAAAAATAGATATTACTATTCTATAATAAGGATGAA".to_string()));
        assert_eq!(reverse_reader.next(), None);
    }

    #[test]
    pub fn parse_anonymous_single_fasta() {
        let filepath = Path::new("./tests/data/anonymous_single_fasta.fa");
        let fasta_reader = Fasta::new(filepath).expect("Error landmarking FASTA file.");
        
        let sequences = fasta_reader.get_sequences();
        assert_eq!(sequences.len(), 1);
        assert_eq!(sequences[0].get_sequence_id(), "Unknown-1");
        assert_eq!(sequences[0].get_sequence_length(), 106);
        
        let mut forward_reader = sequences[0].build_forward_reader().expect("Error building forward reader");

        assert_eq!(forward_reader.next(), Some("AAGTAGGAATAATATCTTATCATTATAGATAAAAACCTTCTGAATTTGCTTAGTGTGTAT".to_string()));
        assert_eq!(forward_reader.next(), Some("ACGACTAGACATATATCAGCTCGCCGATTATTTGGATTATTCCCTG".to_string()));
        assert_eq!(forward_reader.next(), None);

        let mut reverse_reader = sequences[0].build_reverse_reader().expect("Error building reverse reader");

        assert_eq!(reverse_reader.next(), Some("GTCCCTTATTAGGTTTATTAGCCGCTCGACTATATACAGATCAGCA".to_string()));
        assert_eq!(reverse_reader.next(), Some("TATGTGTGATTCGTTTAAGTCTTCCAAAAATAGATATTACTATTCTATAATAAGGATGAA".to_string()));
        assert_eq!(reverse_reader.next(), None);
    }

    #[test]
    pub fn parse_multi_fasta() {
        let filepath = Path::new("./tests/data/multi_fasta.fa");
        let fasta_reader = Fasta::new(filepath).expect("Error landmarking FASTA file.");

        let sequences = fasta_reader.get_sequences();
        assert_eq!(sequences.len(), 3);
        assert_eq!(sequences[0].get_sequence_id(), "sequenceID-001");
        assert_eq!(sequences[1].get_sequence_id(), "sequenceID-002");
        assert_eq!(sequences[2].get_sequence_id(), "sequenceID-003");
        assert_eq!(sequences[0].get_sequence_length(), 106);
        assert_eq!(sequences[1].get_sequence_length(), 154);
        assert_eq!(sequences[2].get_sequence_length(), 76);
        
        let mut forward_reader = sequences[0].build_forward_reader().expect("Error building forward reader");

        assert_eq!(forward_reader.next(), Some("AAGTAGGAATAATATCTTATCATTATAGATAAAAACCTTCTGAATTTGCTTAGTGTGTAT".to_string()));
        assert_eq!(forward_reader.next(), Some("ACGACTAGACATATATCAGCTCGCCGATTATTTGGATTATTCCCTG".to_string()));
        assert_eq!(forward_reader.next(), None);

        let mut reverse_reader = sequences[0].build_reverse_reader().expect("Error building reverse reader");

        assert_eq!(reverse_reader.next(), Some("GTCCCTTATTAGGTTTATTAGCCGCTCGACTATATACAGATCAGCA".to_string()));
        assert_eq!(reverse_reader.next(), Some("TATGTGTGATTCGTTTAAGTCTTCCAAAAATAGATATTACTATTCTATAATAAGGATGAA".to_string()));
        assert_eq!(reverse_reader.next(), None);

        let mut forward_reader = sequences[1].build_forward_reader().expect("Error building forward reader");

        assert_eq!(forward_reader.next(), Some("CAGTAAAGAGTGGATGTAAGAACCGTCCGATCTACCAGATGTGATAGAGGTTGCCAGTAC".to_string()));
        assert_eq!(forward_reader.next(), Some("AAAAATTGCATAATAATTGATTAATCCTTTAATATTGTTTAGAATATATCCGTCAGATAA".to_string()));
        assert_eq!(forward_reader.next(), Some("TCCTAAAAATAACGATATGATGGCGGAAATCGTC".to_string()));
        assert_eq!(forward_reader.next(), None);

        let mut reverse_reader = sequences[1].build_reverse_reader().expect("Error building reverse reader");

        assert_eq!(reverse_reader.next(), Some("CTGCTAAAGGCGGTAGTATAGCAATAAAAATCCT".to_string()));
        assert_eq!(reverse_reader.next(), Some("AATAGACTGCCTATATAAGATTTGTTATAATTTCCTAATTAGTTAATAATACGTTAAAAA".to_string()));
        assert_eq!(reverse_reader.next(), Some("CATGACCGTTGGAGATAGTGTAGACCATCTAGCCTGCCAAGAATGTAGGTGAGAAATGAC".to_string()));
        assert_eq!(reverse_reader.next(), None);

        let mut forward_reader = sequences[2].build_forward_reader().expect("Error building forward reader");

        assert_eq!(forward_reader.next(), Some("CTTCAATTACCCTGCTGACGCGAGATACCTTATGCATCGAAGGTAAAGCGATGAATTTAT".to_string()));
        assert_eq!(forward_reader.next(), Some("CCAAGGTTTTAATTTG".to_string()));
        assert_eq!(forward_reader.next(), None);

        let mut reverse_reader = sequences[2].build_reverse_reader().expect("Error building reverse reader");

        assert_eq!(reverse_reader.next(), Some("GTTTAATTTTGGAACC".to_string()));
        assert_eq!(reverse_reader.next(), Some("TATTTAAGTAGCGAAATGGAAGCTACGTATTCCATAGAGCGCAGTCGTCCCATTAACTTC".to_string()));
        assert_eq!(reverse_reader.next(), None);
    }

    #[test]
    pub fn parse_multi_fasta_extra_linefeeds() {
        let filepath = Path::new("./tests/data/multi_fasta_extra_linefeeds.fa");
        let fasta_reader = Fasta::new(filepath).expect("Error landmarking FASTA file.");

        let sequences = fasta_reader.get_sequences();
        assert_eq!(sequences.len(), 3);
        assert_eq!(sequences[0].get_sequence_id(), "sequenceID-001");
        assert_eq!(sequences[1].get_sequence_id(), "sequenceID-002");
        assert_eq!(sequences[2].get_sequence_id(), "sequenceID-003");
        assert_eq!(sequences[0].get_sequence_length(), 106);
        assert_eq!(sequences[1].get_sequence_length(), 154);
        assert_eq!(sequences[2].get_sequence_length(), 76);
        
        let mut forward_reader = sequences[0].build_forward_reader().expect("Error building forward reader");

        assert_eq!(forward_reader.next(), Some("AAGTAGGAATAATATCTTATCATTATAGATAAAAACCTTCTGAATTTGCTTAGTGTGTAT".to_string()));
        assert_eq!(forward_reader.next(), Some("ACGACTAGACATATATCAGCTCGCCGATTATTTGGATTATTCCCTG".to_string()));
        assert_eq!(forward_reader.next(), None);

        let mut reverse_reader = sequences[0].build_reverse_reader().expect("Error building reverse reader");

        assert_eq!(reverse_reader.next(), Some("GTCCCTTATTAGGTTTATTAGCCGCTCGACTATATACAGATCAGCA".to_string()));
        assert_eq!(reverse_reader.next(), Some("TATGTGTGATTCGTTTAAGTCTTCCAAAAATAGATATTACTATTCTATAATAAGGATGAA".to_string()));
        assert_eq!(reverse_reader.next(), None);

        let mut forward_reader = sequences[1].build_forward_reader().expect("Error building forward reader");

        assert_eq!(forward_reader.next(), Some("CAGTAAAGAGTGGATGTAAGAACCGTCCGATCTACCAGATGTGATAGAGGTTGCCAGTAC".to_string()));
        assert_eq!(forward_reader.next(), Some("AAAAATTGCATAATAATTGATTAATCCTTTAATATTGTTTAGAATATATCCGTCAGATAA".to_string()));
        assert_eq!(forward_reader.next(), Some("TCCTAAAAATAACGATATGATGGCGGAAATCGTC".to_string()));
        assert_eq!(forward_reader.next(), None);

        let mut reverse_reader = sequences[1].build_reverse_reader().expect("Error building reverse reader");

        assert_eq!(reverse_reader.next(), Some("CTGCTAAAGGCGGTAGTATAGCAATAAAAATCCT".to_string()));
        assert_eq!(reverse_reader.next(), Some("AATAGACTGCCTATATAAGATTTGTTATAATTTCCTAATTAGTTAATAATACGTTAAAAA".to_string()));
        assert_eq!(reverse_reader.next(), Some("CATGACCGTTGGAGATAGTGTAGACCATCTAGCCTGCCAAGAATGTAGGTGAGAAATGAC".to_string()));
        assert_eq!(reverse_reader.next(), None);

        let mut forward_reader = sequences[2].build_forward_reader().expect("Error building forward reader");

        assert_eq!(forward_reader.next(), Some("CTTCAATTACCCTGCTGACGCGAGATACCTTATGCATCGAAGGTAAAGCGATGAATTTAT".to_string()));
        assert_eq!(forward_reader.next(), Some("CCAAGGTTTTAATTTG".to_string()));
        assert_eq!(forward_reader.next(), None);

        let mut reverse_reader = sequences[2].build_reverse_reader().expect("Error building reverse reader");

        assert_eq!(reverse_reader.next(), Some("GTTTAATTTTGGAACC".to_string()));
        assert_eq!(reverse_reader.next(), Some("TATTTAAGTAGCGAAATGGAAGCTACGTATTCCATAGAGCGCAGTCGTCCCATTAACTTC".to_string()));
        assert_eq!(reverse_reader.next(), None);
    }

    #[test]
    pub fn parse_anonymous_multi_fasta() {
        let filepath = Path::new("./tests/data/anonymous_multi_fasta.fa");
        let fasta_reader = Fasta::new(filepath).expect("Error landmarking FASTA file.");

        let sequences = fasta_reader.get_sequences();
        assert_eq!(sequences.len(), 3);
        assert_eq!(sequences[0].get_sequence_id(), "Unknown-1");
        assert_eq!(sequences[1].get_sequence_id(), "Unknown-2");
        assert_eq!(sequences[2].get_sequence_id(), "Unknown-3");
        assert_eq!(sequences[0].get_sequence_length(), 106);
        assert_eq!(sequences[1].get_sequence_length(), 154);
        assert_eq!(sequences[2].get_sequence_length(), 76);
        
        let mut forward_reader = sequences[0].build_forward_reader().expect("Error building forward reader");

        assert_eq!(forward_reader.next(), Some("AAGTAGGAATAATATCTTATCATTATAGATAAAAACCTTCTGAATTTGCTTAGTGTGTAT".to_string()));
        assert_eq!(forward_reader.next(), Some("ACGACTAGACATATATCAGCTCGCCGATTATTTGGATTATTCCCTG".to_string()));
        assert_eq!(forward_reader.next(), None);

        let mut reverse_reader = sequences[0].build_reverse_reader().expect("Error building reverse reader");

        assert_eq!(reverse_reader.next(), Some("GTCCCTTATTAGGTTTATTAGCCGCTCGACTATATACAGATCAGCA".to_string()));
        assert_eq!(reverse_reader.next(), Some("TATGTGTGATTCGTTTAAGTCTTCCAAAAATAGATATTACTATTCTATAATAAGGATGAA".to_string()));
        assert_eq!(reverse_reader.next(), None);

        let mut forward_reader = sequences[1].build_forward_reader().expect("Error building forward reader");

        assert_eq!(forward_reader.next(), Some("CAGTAAAGAGTGGATGTAAGAACCGTCCGATCTACCAGATGTGATAGAGGTTGCCAGTAC".to_string()));
        assert_eq!(forward_reader.next(), Some("AAAAATTGCATAATAATTGATTAATCCTTTAATATTGTTTAGAATATATCCGTCAGATAA".to_string()));
        assert_eq!(forward_reader.next(), Some("TCCTAAAAATAACGATATGATGGCGGAAATCGTC".to_string()));
        assert_eq!(forward_reader.next(), None);

        let mut reverse_reader = sequences[1].build_reverse_reader().expect("Error building reverse reader");

        assert_eq!(reverse_reader.next(), Some("CTGCTAAAGGCGGTAGTATAGCAATAAAAATCCT".to_string()));
        assert_eq!(reverse_reader.next(), Some("AATAGACTGCCTATATAAGATTTGTTATAATTTCCTAATTAGTTAATAATACGTTAAAAA".to_string()));
        assert_eq!(reverse_reader.next(), Some("CATGACCGTTGGAGATAGTGTAGACCATCTAGCCTGCCAAGAATGTAGGTGAGAAATGAC".to_string()));
        assert_eq!(reverse_reader.next(), None);

        let mut forward_reader = sequences[2].build_forward_reader().expect("Error building forward reader");

        assert_eq!(forward_reader.next(), Some("CTTCAATTACCCTGCTGACGCGAGATACCTTATGCATCGAAGGTAAAGCGATGAATTTAT".to_string()));
        assert_eq!(forward_reader.next(), Some("CCAAGGTTTTAATTTG".to_string()));
        assert_eq!(forward_reader.next(), None);

        let mut reverse_reader = sequences[2].build_reverse_reader().expect("Error building reverse reader");

        assert_eq!(reverse_reader.next(), Some("GTTTAATTTTGGAACC".to_string()));
        assert_eq!(reverse_reader.next(), Some("TATTTAAGTAGCGAAATGGAAGCTACGTATTCCATAGAGCGCAGTCGTCCCATTAACTTC".to_string()));
        assert_eq!(reverse_reader.next(), None);
    }
}
