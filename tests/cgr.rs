use chaos_game_representation::{
    BufferedChaosGameRepresentation, 
    ChaosGameRepresentation, 
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
pub fn cgr_test() {
    let filename = "./data/egfr/egfr.fa";
    //let filename = "./data/grch37/hg19.fa";

    let filepath = Path::new(filename);
    let outdir = Path::new("./data/figures");
    let cgr = ChaosGameRepresentation::from_fasta_file(filepath);
    cgr.plot(outdir, 1024, 30).expect("Error in plotting");
}


#[test]
pub fn buffered_cgr_test() {
    //let filepath = Path::new("./tests/data/multi_fasta.fa");
    let filepath = Path::new("./tests/data/single_fasta.fa");

    let cgr_buf = BufferedChaosGameRepresentation::new(filepath);
    //let filepath_hdf5 = Path::new("./tests/data/multi.h5");
    let filepath_hdf5 = Path::new("./tests/data/single.h5");
    cgr_buf.build_cgrs_and_write_to_hdf5(filepath_hdf5, 10).expect("Error in HDF5 file creation.");
}