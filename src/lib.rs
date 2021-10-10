extern crate pyo3;
extern crate bio;
extern crate flate2;
use std::string::String;
use std::collections::HashMap;
use std::convert::TryInto;
use bio::io::fastq;
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use pyo3::exceptions::PyValueError;
mod utils;
mod bed;
mod bed12;



#[pyfunction]
/// ---
/// 
/// Function to calculate how many bases and how many record in the fastq file
/// 
/// Args:
///     filename: str
///         full path to the fastq file, can be gz zipped file
/// 
/// return:
///     read_count: int
///         how many reads are in the fastq file?
///     base_count: int
///         how many bases in total are in the fastq file?
fn fq_stat(filename: String) -> PyResult<(usize, usize)>{
    let reader = fastq::Reader::new(utils::read_file(&filename));
    let mut basecount = 0;
    let mut readcount = 0;
    for result in reader.records(){
        let record = result.expect("Error during fastq record parsing");
        let seq = record.seq();
        let len = seq.len();
        basecount += len;
        readcount += 1;
    }
    Ok((readcount, basecount))
}


#[pyfunction]
#[text_signature = "(seq, k, /)"]
/// --
///
/// counting kmer from the input sequence
/// 
/// Args:
///     seq: str
///         read sequence
///     k: int
///         kmer size
/// 
/// Return:
///     kmer_dict: dict
///         with kmer as key and, kmer count as value
fn kmer_counter(seq: String, k: usize) -> PyResult<HashMap<String, usize>>{
    let mut kmer_count: HashMap<String, usize> = HashMap::new();
    let seq_len: usize = seq.len().try_into().unwrap();
    if (seq_len < k) {
        Err(PyValueError::new_err("k is smaller than sequence length"))
    } else {
		for i in 0..(seq_len - k + 1){
			//let kmer = utils::substring(seq, i, i+k);
			let kmer = seq[i..(i+k)].to_string();
			*kmer_count.entry(kmer).or_insert(0) += 1;
		}
		Ok(kmer_count)
	}
}



/// A Python module implemented in Rust.
#[pymodule]
fn biotools_lib(py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(fq_stat))?;
    m.add_wrapped(wrap_pyfunction!(kmer_counter))?;
    m.add_class::<bed::BedRecord>()?;
    m.add_class::<bed12::Bed12Record>()?;
    Ok(())
}
