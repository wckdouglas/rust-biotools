use std::string::String;
use std::format;
use pyo3::prelude::*;
use pyo3::{Py, PyResult};
use std::cmp::{max,min};

#[pyclass]
pub struct BedRecord {
    #[pyo3(get, set)]
    pub chrom: String,
    #[pyo3(get, set)]
    pub start: i32,
    #[pyo3(get, set)]
    pub end: i32,
    #[pyo3(get, set)]
    pub name: String,
    #[pyo3(get, set)]
    pub strand: String
}


#[pymethods]
impl BedRecord {
    #[new]
    pub fn new(bedline: String) -> Self {
        println!("{}",bedline);
        let v: Vec<&str> = bedline.trim().split("\t").collect();
        let mut it = v.iter();
        let chrom = it.next().expect("no chrom field").to_string();
        let start: i32 = it.next().expect("no start field").parse().unwrap();
        let end: i32 = it.next().expect("no end field").parse().unwrap();
        let name = it.next().expect("no name field").to_string();
        let _ = it.next().expect("no score field").to_string();
        let strand = it.next().expect("no strand field").to_string();
        Self{
            chrom: chrom,
            start: start,
            end: end,
            name: name,
            strand: strand
        }
    }

    #[getter]
    fn coordinate(&self) -> PyResult<String> {
        Ok(format!("{}:{}-{}", self.chrom, self.start, self.end))
    }

    pub fn overlap(&self, start: i32, end: i32) -> PyResult<bool>{
        if max(self.start, start) <= min(self.end, end) {
            Ok(true)
        }else{
            Ok(false)
        }
    }
}