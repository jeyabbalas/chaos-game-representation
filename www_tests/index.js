// Serve the HDF5 file using $> http-server -p 8000 --cors
import * as hdf5 from "hdf5";

await hdf5.ready;

let response = await fetch("http://127.0.0.1:8000/hg19_cgr.h5");
let ab = await response.arrayBuffer();
hdf5.FS.writeFile("hg19_cgr.h5", new Uint8Array(ab));
let f = new hdf5.File("hg19_cgr.h5", "r");
            
console.log(f.keys());

let chr21 = f.get("chr21");
console.log(chr21.keys());

let chr21_seq = chr21.get("sequence");
let chr21_fwd = chr21.get("forward_cgr");
let chr21_bwd = chr21.get("backward_cgr");

console.log(chr21_seq.shape);

let idx = 16482558
console.log("Reference allele: ");
console.log(chr21_seq.slice( [[idx-1,idx]] ));
console.log(chr21_seq.slice( [[idx-2,idx+1]] ));
console.log("\nForward coordinate: ");
console.log(chr21_fwd.slice( [[idx-2,idx-1], []] ));
console.log("\nBackward coordinate: ");
console.log(chr21_bwd.slice( [[idx,idx+1], []] ));

f.close();
hdf5.Module.FS.unlink("hg19_cgr.h5");
