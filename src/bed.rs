use std::string::String;
use std::format;
use pyo3::prelude::*;
use pyo3::PyResult;

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
        let mut it = bedline.split("\t");
        let chrom = it.next().expect("no chrom field").to_string();
        let start: i32 = it.next().expect("no start field").parse().unwrap();
        let end: i32 = it.next().expect("no end field").parse().unwrap();
        let name = it.next().expect("no name field").to_string();
        let _ = it.next().expect("no score field").to_string();
        let strand = it.next().expect("no strand field").to_string();
        BedRecord{chrom: chrom, 
                start: start, 
                end: end,
                name: name,
                strand: strand}
    }

    #[getter]
    fn coordinate(&self) -> PyResult<String> {
        Ok(format!("{}:{}-{}", self.chrom, self.start, self.end))
    }

}