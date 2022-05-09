use chaos_game_representation::{
    BufferedChaosGameRepresentation, 
    Point
};

use std::path::Path;


const EPSILON: f64 = 0.0000001;

#[test]
pub fn add_test() {
    let sum = Point { x: 0.5, y: 1.3 } + Point { x: 2.0, y: 3.7 };
    assert!(((sum.x - 2.5) as f64).abs() < EPSILON);
    assert!(((sum.y - 5.0) as f64).abs() < EPSILON);
}


#[test]
pub fn subtract_test() {
    let sub = Point { x: 0.5, y: 3.7 } - Point { x: 2.0, y: 1.3 };
    assert!(((sub.x + 1.5) as f64).abs() < EPSILON);
    assert!(((sub.y - 2.4) as f64).abs() < EPSILON);
}


#[test]
pub fn divide_test() {
    let div = Point { x: 2.0, y: 0.5 } / 2.0;
    assert!(((div.x - 1.0) as f64).abs() < EPSILON);
    assert!(((div.y - 0.25) as f64).abs() < EPSILON);
}


#[test]
pub fn buffered_cgr_test() {
    let filepath = Path::new("./tests/data/multi_fasta.fa");
    //let filepath = Path::new("./tests/data/single_fasta.fa");

    let cgr_buf = BufferedChaosGameRepresentation::new(filepath);
    let output_dir = Path::new("./tests/data/output/");
    cgr_buf.build_cgrs_and_write_to_binary_file(output_dir).expect("Error in CGR creation.");
}


#[test]
#[should_panic]
pub fn test_invalid_sequence_id() {
    let filepath = Path::new("./tests/data/multi_fasta.fa");

    let cgr_buf = BufferedChaosGameRepresentation::new(filepath);
    let sequence_ids = vec![
        "Boo!"
    ];
    let output_dir = Path::new("./tests/data/output/");
    cgr_buf.build_select_cgrs_and_write_to_binary_file(output_dir, sequence_ids).expect("Error in CGR creation.");
}


#[test]
pub fn test_temp() {
    use std::fs::File;
    use std::io::prelude::*;

    let mut file = File::create("test_float").expect("Error creating file");
    let floats = [0.252627_f64, 0.599934, 0.999999, 0.01010101];
    let mut msg = Vec::<u8>::with_capacity(floats.len() * 8);
    for float in floats.iter() {
        msg.extend_from_slice(&float.to_be_bytes());
    }
    file.write_all(&msg).expect("Trouble writing");
}
