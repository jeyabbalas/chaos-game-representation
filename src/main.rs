use std::{env, process};
use chaos_game_representation::Config;


fn main() {
    let config = Config::new(env::args()).unwrap_or_else(|err| {
        eprintln!("Error while parsing input arguments: {}", err);
        process::exit(1);
    });

    if let Err(e) = chaos_game_representation::run(config) {
        eprintln!("Error while reading input file: {}.", e);
        process::exit(1);
    }
}
