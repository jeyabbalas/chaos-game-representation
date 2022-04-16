// Serve the HDF5 file using $> http-server -p 8000 --cors
import * as hdf5 from "hdf5";

await hdf5.ready;

let response = await fetch("http://localhost:8000/hg19_cgr.h5");
let ab = await response.arrayBuffer();
hdf5.FS.writeFile("hg19_cgr.h5", new Uint8Array(ab));
let f = new hdf5.File("hg19_cgr.h5", "r");
            
console.log(f.keys());

let chr21 = f.get("chr21");
console.log(chr21.keys());

let chr21_seq = chr21.get("sequence");
console.log(chr21_seq.shape);
console.log(chr21_seq.slice([[18002381,18002387]]));

f.close();
hdf5.Module.FS.unlink("hg19_cgr.h5");
