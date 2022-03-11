use std::path::Path;


fn main() {
    let fasta_filepath: &Path = Path::new("./data/egfr/egfr.fa.gz");
    //let fasta_filepath: &Path = Path::new("./data/grch37/hg19.fa.gz");
    chaos_game_representation::parse_fasta_file(fasta_filepath);
    //chaos_game_representation::temp();
}