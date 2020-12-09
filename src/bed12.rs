use std::string::String;
use std::format;
use pyo3::prelude::*;
use pyo3::{Py, PyResult};
use std::cmp::{max,min};

#[pyclass]
pub struct Bed12Record {
    #[pyo3(get, set)]
    pub chrom: String,
    #[pyo3(get, set)]
    pub start: i32,
    #[pyo3(get, set)]
    pub end: i32,
    #[pyo3(get, set)]
    pub name: String,
    #[pyo3(get, set)]
    pub score: i32,
    #[pyo3(get, set)]
    pub strand: String,
    #[pyo3(get, set)]
    pub coding_start: i32,
    #[pyo3(get, set)]
    pub coding_end: i32,
    #[pyo3(get, set)]
    pub block_count: i32,
    #[pyo3(get, set)]
    pub block_starts: Vec<i32>,
    #[pyo3(get, set)]
    pub block_ends: Vec<i32>,
}


#[pymethods]
impl Bed12Record {
    #[new]
    pub fn new(bedline: String) -> Self {
        println!("{}",bedline);
        let v: Vec<&str> = bedline.trim().split("\t").collect();
        let mut it = v.iter();
        let chrom = it.next().expect("no chrom field").to_string();
        let start: i32 = it.next().expect("no start field").parse().unwrap();
        let end: i32 = it.next().expect("no end field").parse().unwrap();
        let name = it.next().expect("no name field").to_string();
        let score = it.next().expect("no score field").parse().unwrap();
        let strand = it.next().expect("no strand field").to_string();
        let coding_start: i32 = it.next().expect("no thickStart field").parse().unwrap();
        let coding_end: i32 = it.next().expect("no thickEnd field").parse().unwrap();
        let _ = it.next().expect("no thickEnd field").to_string();
        let block_count: i32 = it.next().expect("no BlockCount field").parse().unwrap();
        let block_sizes: Vec<i32> = it.next().expect("no blockSizes field").to_string().trim_end_matches(',').split(",").map(|x| x.parse::<i32>().unwrap()).collect();
        let mut block_starts: Vec<i32> = it.next().expect("no blockStarts field").to_string().trim_end_matches(',').split(",").map(|x| x.parse::<i32>().unwrap()).collect();
        block_starts = block_starts.into_iter().map(|x| x + start).collect();
        let (mut block_starts, mut block_ends) = adding(block_starts, block_sizes);
        if strand == "-" {
            block_starts.reverse();
            block_ends.reverse();
        }
        Self{
            chrom: chrom,
            start: start,
            end: end,
            name: name,
            score: score,
            strand: strand,
            coding_start: coding_start,
            coding_end: coding_end,
            block_count: block_count,
            block_ends: block_ends,
            block_starts: block_starts,
        }
    }


    #[getter]
    fn coordinate(&self) -> PyResult<String> {
        Ok(format!("{}:{}-{}", self.chrom, self.start, self.end))
    }

    /// ---
    /// 
    /// Function to check overlap between ranges
    /// 
    /// Args:
    ///     start (int): leftmost position of the range
    ///     end (int): rightmost position of the range
    /// 
    /// return:
    ///     boolean: True if overlap
    pub fn overlap(&self, start: i32, end: i32) -> PyResult<bool>{
        for (s, e) in self.block_starts.iter().zip(self.block_ends.iter()){
            if max(s, &start) <= min(e, &end) {
                return Ok(true);
            }
        }
        return Ok(false);
    }

}

fn adding(starts: Vec<i32>, sizes: Vec<i32>) -> (Vec<i32>,Vec<i32>){
    let mut res = vec![];
    for (start, size) in starts.iter().zip(sizes.iter()){
        let sum = start + size;
        res.push(sum);
    }
    (starts, res)
}