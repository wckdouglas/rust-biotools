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
    pub exons: Vec<Exon>,
}

#[pyclass]
pub struct Gene {
    pub chrom: String,
    pub strand: String,
    pub transcripts: HashMap<String, &Transcript>
}

impl Gene {
    pub fn add_transcript(&mut self, rnaId: &String, transcript: &Transcript){
        self.transcripts.insert((&rnaId).to_string(), 
                                transcript);
    }
}

pub struct Transcriptome{
    transcriptome: HashMap<String, Gene>,
}


impl Transcriptome {
    pub fn new(refflat: String) -> Self {
        let mut transcriptome: HashMap<String, Gene> = HashMap::new();

        if let Ok(lines) = utils::read_lines(refflat){
            for line in lines {
                if let l = line.unwrap().to_string() {
                    let fields: Vec<&str> = l.trim().split('\t').collect();
                    let (gene_name, rnaId, chrom, strand, transcript) = Self::parse_transcript(fields);
                    let mut transcripts: HashMap<String, Transcript> = HashMap::new();
                    transcripts.insert(rnaId, transcript);
                    let mut gene = Gene{chrom: chrom,
                                        strand: strand,
                                        transcripts: transcripts};
                    if ! transcriptome.contains_key(&gene_name){
                        transcriptome.insert(gene_name, gene);
                    } else {
                        transcriptome
                            .entry(gene_name)
                            .or_insert(gene)
                            .add_transcript( &rnaId, transcript ) ;
                    }
                }
            }
        }
        Self{transcriptome: transcriptome}
    }

    fn parse_transcript(fields: Vec<&str>) -> (String, String, String, String, Transcript) {
        let gene_name = fields[0].to_string();
        let rnaId = fields[1].to_string();
        let chrom = fields[2].to_string();
        let strand = fields[3].to_string();
        let strand = fields[3].to_string();
        let tx_start: i32 = fields[4].parse().unwrap();
        let tx_end: i32 = fields[5].parse().unwrap();
        let cds_start: i32 = fields[6].parse().unwrap();
        let cds_end: i32 = fields[7].parse().unwrap();
        let exon_count: usize = fields[8].parse().unwrap();
        let exon_starts = fields[9].trim_end_matches(',').to_string();
        let exon_ends = fields[10].trim_end_matches(',').to_string();
        let exon_vec = Self::parse_exons(exon_starts, exon_ends, cds_start, cds_end, &strand);
        assert_eq!(exon_count, exon_vec.len(), "Not all exons being parsed");
        let transcript = Transcript{ start: tx_start,
                                    end: tx_end,
                                    exons: exon_vec};
        (gene_name, rnaId, chrom, strand, transcript)
    }

    fn parse_exons(exon_starts: String, exon_ends: String, cds_start: i32, cds_end: i32, strand: &String) -> Vec<Exon>{
        /*

        Parsing the last two columns from refFlat file to make exons:

        Forward strand transcripts:

              CDS                                                CDE
               |->                                               ||
        5' |=========|---|========|---------|========|--------|=====|> 3'
            exon 1          exon 2            exon 3           exon 4

              CDE                                               CDS
               ||                                               <-|
        3' <|=========|---|========|---------|========|--------|=====| 5'
            exon 4          exon 3            exon 2           exon 1

        */
        let mut exon_vec = vec![];
        let mut exon_start_vec: Vec<&str> = exon_starts.split(',').collect();
        let mut exon_end_vec: Vec<&str> = exon_ends.split(',').collect();
        assert_eq!(exon_end_vec.len(), exon_start_vec.len());
        let mut is_coding = true;
        if cds_start == cds_end{
            is_coding = false;
        }

        if strand == "-" {
            // flip the exon order if reverse strand transcript
            exon_end_vec.reverse();
            exon_start_vec.reverse();
        }

        let zipped = exon_start_vec.iter().zip(exon_end_vec.iter());

        let mut coding_exon = false; // is a coding exon? if it is after cds, or contain cds, or contain cde
        let mut contain_cds = false; // containing coding start site?
        let mut contain_cde = false; // containing coding end site?
        for (es, ee) in zipped {
            let exon_start = es.parse().unwrap();
            let exon_end = ee.parse().unwrap();
            
            if is_coding {
                if exon_start <= cds_start && cds_start <= exon_end {
                    contain_cds = true;
                    coding_exon = true;
                    if exon_start <= cds_end && cds_end <= exon_end {
                        /* contain CDS and CDE 
                        CDS           CDE
                          |->         ||
                        |===================|-----------|============|
                            (this exon)
                        */
                        contain_cde = true;
                        let exon = Exon::new(exon_start,
                                            exon_end,
                                            coding_exon,
                                            contain_cds,
                                            contain_cde);
                        exon_vec.push(exon);
                        coding_exon = false; //no more coding exon after this
                    } else{
                        /* contain CDS but no 
                        CDS                                     CDE
                          |->                                    ||
                        |===================|--------------|===========|
                            (this exno)
                        */
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
                        /* doesn't contain CDS but contain CDE
                        CDS                                 CDE
                          |->                                ||
                        |===================|-----------|============|
                                                           (this exon)
                        */
                        contain_cde = true;
                        let exon = Exon::new(exon_start,
                                            exon_end,
                                            coding_exon,
                                            contain_cds,
                                            contain_cde);
                        exon_vec.push(exon);
                        coding_exon = false;
                    }else{
                        /* doesn't contain CDS but contain CDE
                        CDS                                                     CDE
                          |->                                                   ||
                        |===================|-----------|============|-------|========|
                                                           (this exon)
                        */
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
                // non coding RNA
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
