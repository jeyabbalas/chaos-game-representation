use chaos_game_representation::BufferedChaosGameRepresentation;

use std::path::Path;


fn main() {
    let filepath = Path::new("./data/grch37/hg19.fa");
    let cgr_buf = BufferedChaosGameRepresentation::new(filepath);
    let select = vec![
        "chr21"
    ];
    //let filepath_hdf5 = Path::new("./data/grch37/hg19_cgr.h5");
    //cgr_buf.build_cgrs_and_write_to_hdf5(filepath_hdf5, 1_000_000).expect("Error in HDF5 file creation.");
}