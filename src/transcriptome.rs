use crate::utils;
use std::string::String;
use std::format;
use std::vec::Vec;
use pyo3::prelude::*;
use pyo3::{Py, PyResult};
use std::cmp::{max,min};

pub struct Exon {
    pub start: i32,
    pub end: i32,
    pub is_coding: bool,
    pub contain_cds: bool,
    pub contain_cde: bool
}

pub struct Transcript {
    pub chrom: String,
    pub start: i32,
    pub end: i32,
    pub strand: String,
    pub exons: Vec<Exon>
}

pub struct Transcriptome{
    Transcripts: String
}

impl Transcriptome {
    pub fn new(refflat: String) -> () {
        if let Ok(lines) = utils::read_lines(refflat){
            for line in lines {
                if let Ok(ip) = line {
                    println!("{}", ip);
                }
            }
        }
    }

}
