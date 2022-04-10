use hdf5::{self, types::FixedAscii};
use ndarray::{Array1, Array2};


fn main() {
    let file = hdf5::File::create("chunked.h5")
        .expect("Error creating HDF5 file.");
    
    let (nx, ny) = (100, 2);
    let arr = Array2::from_shape_fn((nx, ny), |(_, _)| 1.0);

    let ds = file.new_dataset::<f64>()
        .chunk((1, nx, ny))
        .shape((19, nx, ny))
        .deflate(3)
        .create("cgr")
        .expect("Error creating HDF5 dataset.");
    
    ds.write_slice(&arr, (0, .., ..))
        .expect("Error writing to HDF5");
    
    // ds.resize((10, nx, ny)).expect("Something");
    ds.write_slice(&arr, (2, .., ..))
        .expect("Error writing to HDF5");

    
    let sequence = "GATTACA";
    let sequence_bytes: Vec<FixedAscii<1>> = sequence.to_ascii_uppercase()
        .chars()
        .map(|base| FixedAscii::<1>::from_ascii(&[base as u8]).unwrap())
        .collect();
    let carr = Array1::from(sequence_bytes);

    let ds_seq = file.new_dataset::<FixedAscii<1>>()
        .chunk(7)
        .shape(21)
        .deflate(9)
        .create("seq")
        .expect("Char array dataset creation error.");
    
    ds_seq.write_slice(&carr, 0..7).expect("Error writing char array");
    ds_seq.write_slice(&carr, 14..21).expect("Error writing char array");
    
}