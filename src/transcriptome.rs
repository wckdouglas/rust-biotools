use crate::utils;
use std::string::String;
use std::vec::Vec;
use std::collections::HashMap;
use pyo3::prelude::*;
use pyo3::{Py, PyResult};

#[pyclass]
pub struct Exon {
    pub start: i32,
    pub end: i32,
    pub exon_length: i32,
    pub coding_exon: bool,
    pub contain_cds: bool,
    pub contain_cde: bool,
}

impl Exon {
    fn new(start: i32, end: i32, 
            coding_exon: bool, 
            contain_cds: bool, 
            contain_cde: bool) -> Self{
        let exon_length =  end - start;
        Self{start: start,
            end: end,
            exon_length: exon_length,
            coding_exon: coding_exon,
            contain_cde: contain_cde,
            contain_cds: contain_cds}
    }
}

#[pyclass]
pub struct Transcript {
    pub start: i32,
    pub end: i32,
    pub exons: Vec<Py<Exon>>,
}

#[pyclass]
pub struct Gene {
    pub chrom: String,
    pub strand: String,
    pub transcripts: HashMap<String, Py<Transcript>>
}

pub struct Transcriptome{
    genes: HashMap<String, Py<Gene>>,
}


impl Transcriptome {
    pub fn new(refflat: String) -> () {
        let mut transcriptome: HashMap<String, Py<Transcript>> = HashMap::new();

        if let Ok(lines) = utils::read_lines(refflat){
            for line in lines {
                if let l = line.unwrap().to_string() {
                    let fields: Vec<&str> = l.trim().split('\t').collect();
                    Self::parse_transcript(fields);
                }
            }
        }
    }

    fn parse_transcript(fields: Vec<&str>) -> (){
        let gene_name = fields[0].to_string();
        let rnaId = fields[1].to_string();
        let chrom = fields[2].to_string();
        let strand = fields[3].to_string();
        let tx_start: i32 = fields[4].parse().unwrap();
        let tx_end: i32 = fields[5].parse().unwrap();
        let cds_start: i32 = fields[6].parse().unwrap();
        let cds_end: i32 = fields[7].parse().unwrap();
        let exon_count: u32 = fields[8].parse().unwrap();
        let exon_starts = fields[8].trim_end_matches(',').to_string();
        let exon_ends = fields[9].trim_end_matches(',').to_string();
        let exon_vec = Self::parse_exons(exon_starts, exon_ends, cds_start, cds_end, strand);
        println!("{}", exon_vec[0].end);

        /*
        if 
        *transcriptome.entry(gene_name) = ;
        */
    }

    fn parse_exons(exon_starts: String, exon_ends: String, cds_start: i32, cds_end: i32, strand: String) -> Vec<Exon>{
        let mut exon_vec = vec![];
        let mut exon_start_vec: Vec<&str> = exon_starts.split(',').collect();
        let mut exon_end_vec: Vec<&str> = exon_ends.split(',').collect();
        let mut is_coding = true;
        if cds_start == cds_end{
            is_coding = false;
        }

        if strand == "-" {
            exon_end_vec.reverse();
            exon_start_vec.reverse();
        }

        let zipped = exon_start_vec.iter().zip(exon_end_vec.iter());

        let mut coding_exon = false;
        let mut contain_cds = false;
        let mut contain_cde = false;
        for (es, ee) in zipped {
            let exon_start = es.parse().unwrap();
            let exon_end = ee.parse().unwrap();
            
            if is_coding {
                if exon_start <= cds_start && cds_start <= exon_end {
                    contain_cds = true;
                    coding_exon = true;
                    if exon_start <= cds_end && cds_end <= exon_end {
                        let contain_cde = true;
                        let exon = Exon::new(exon_start,
                                            exon_end,
                                            coding_exon,
                                            contain_cds,
                                            contain_cde);
                        exon_vec.push(exon);
                        coding_exon = false;
                    } else{
                        contain_cde = false;
                        let exon = Exon::new(exon_start,
                                            exon_end,
                                            coding_exon,
                                            contain_cds,
                                            contain_cde);
                        exon_vec.push(exon);
                    }
                }else {
                    contain_cds = false;
                    if exon_start <= cds_end && cds_end <= exon_end{
                        contain_cde = true;
                        let exon = Exon::new(exon_start,
                                            exon_end,
                                            coding_exon,
                                            contain_cds,
                                            contain_cde);
                        exon_vec.push(exon);
                        coding_exon = false;
                    }else{
                        contain_cde = false;
                        let exon = Exon::new(exon_start,
                                            exon_end,
                                            coding_exon,
                                            contain_cds,
                                            contain_cde);
                        exon_vec.push(exon);
                    }
                }

            } else{
                contain_cds = false;
                contain_cde = false;
                coding_exon = false;
                let exon = Exon::new(exon_start,
                                    exon_end,
                                    coding_exon,
                                    contain_cds,
                                    contain_cde);
                exon_vec.push(exon);
            }
        }
        exon_vec
    }
}
