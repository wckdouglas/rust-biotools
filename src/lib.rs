extern crate bio;
extern crate flate2;
use std::string::String;
use std::io::BufReader;
use std::fs;
use bio::io::fastq;
use bio::io::fastq::FastqRead;
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use flate2::bufread;

fn get_fastq_reader(path: &String) -> Box<::std::io::Read> {
    // borrowed from 
    // https://github.com/sndrtj/fastq-count/blob/master/src/main.rs
    if path.ends_with(".gz") {
        let f = fs::File::open(path).unwrap();
        Box::new(bufread::MultiGzDecoder::new(BufReader::new(f)))
    } else {
        Box::new(fs::File::open(path).unwrap())
    }
}

#[pyfunction]
fn readfq(filename: String) -> PyResult<(usize, usize)>{
	let reader = fastq::Reader::new(get_fastq_reader(&filename));
    let mut basecount = 0;
    let mut readcount = 0;
	for record in reader.records(){
        let len = record.unwrap().seq().len();
        basecount += len;
        readcount += 1;
    }
	Ok((readcount, basecount))
}


/// A Python module implemented in Rust.
#[pymodule]
fn biotools(py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(readfq, m)?)?;
    Ok(())
}
