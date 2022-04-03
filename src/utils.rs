use flate2::bufread;
use std::fs;
use std::io::prelude::*;
use std::io::{BufReader, Lines, Result};
use std::path::Path;
use std::string::String;

pub fn read_file(path: &String) -> Box<::std::io::Read> {
    // borrowed from
    // https://github.com/sndrtj/fastq-count/blob/master/src/main.rs
    if path.ends_with(".gz") {
        let f = fs::File::open(path).unwrap();
        Box::new(bufread::MultiGzDecoder::new(BufReader::new(f)))
    } else {
        Box::new(fs::File::open(path).unwrap())
    }
}

// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
pub fn read_lines<P>(filename: P) -> Result<Lines<BufReader<fs::File>>>
where
    P: AsRef<Path>,
{
    let file = fs::File::open(filename)?;
    Ok(BufReader::new(file).lines())
}
