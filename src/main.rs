use chaos_game_representation::BufferedChaosGameRepresentation;

use std::path::Path;


fn main() {
    //let filepath = Path::new("./data/grch37/hg19.fa");
    let filepath = Path::new("./tests/data/single_fasta.fa");
    let cgr_buf = BufferedChaosGameRepresentation::new(filepath);
    let sequence_ids = vec![
        "chr21"
    ];
    //let filepath_hdf5 = Path::new("./data/grch37/hg19_cgr.h5");
    let filepath_hdf5 = Path::new("./tests/data/single.h5");
    cgr_buf.build_select_cgrs_and_write_to_hdf5(filepath_hdf5, 1_000_000, sequence_ids)
        .expect("Error building CGRs and writing it to HDF5.");
}