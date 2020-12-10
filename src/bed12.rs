use std::string::String;
use std::format;
use pyo3::prelude::*;
use pyo3::{Py, PyResult};
use std::cmp::{max,min};

#[pyclass]
pub struct Exon {
    #[pyo3(get, set)]
    pub tstart: i32, // transcript start
    #[pyo3(get, set)]
    pub tend: i32, // transcript end
    #[pyo3(get, set)]
    pub gstart: i32, // genome start
    #[pyo3(get, set)]
    pub gend: i32, // genome end
    #[pyo3(get, set)]
    pub exon_size: i32, // size of exon
}


impl Exon {
    fn new(gstart: i32, gend: i32, tstart: i32) -> Self {
        let exon_size = gend - gstart;
        let tend = tstart + exon_size;
        Self{tstart: tstart, 
            tend: tend,
            gstart: gstart,
            gend: gend,
            exon_size: exon_size}
    }

    fn gcontains(self, gpos: i32) -> bool {
        /// contain genomic coordinate
        if exon.gstart <= gpos && gpos <= self.end {
            return true
        }
        else {
            return false
        }
    }

    fn tcontains(self, tpos: i32) -> bool {
        /// contain transcript coordinate
        if exon.tstart <= tpos && tpos <= self.tend {
            return true
        } else {
            return false
        }
    }
}


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
    #[pyo3(get)]
    pub exons: Vec<Exon>, //exon 1 always at Vec[0], reversed for reverse strand
    #[pyo3(get, set)]
    pub transcript_length: i32,
}


#[pymethods]
impl Bed12Record {
    #[new]
    pub fn new(bedline: String) -> Self {
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
        let transcript_length: i32 = block_sizes.iter().sum();
        let mut block_starts: Vec<i32> = it.next().expect("no blockStarts field").to_string().trim_end_matches(',').split(",").map(|x| x.parse::<i32>().unwrap()).collect();
        block_starts = block_starts.into_iter().map(|x| x + start).collect();
        let (mut block_starts, mut block_ends) = adding(block_starts, block_sizes);
        if strand == "-" {
            block_starts.reverse();
            block_ends.reverse();
        }
        let exon_vec = parse_exons(block_starts, block_ends);
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
            exons: exon_vec,
            transcript_length: transcript_length,
        }
    }


    fn parse_exons(exon_start_vec: Vec<i32>, exon_end_vec: Vec<i32>) -> Vec<Exon>{
        /*

        Parsing the last two columns from bed12 file to make exons:

        Forward strand transcripts:
        5' |=========|---|========|---------|========|--------|=====|> 3'
            exon 1          exon 2            exon 3           exon 4

        3' <|=========|---|========|---------|========|--------|=====| 5'
            exon 4          exon 3            exon 2           exon 1

        */
        let mut exon_vec = vec![];
        assert_eq!(exon_end_vec.len(), exon_start_vec.len());
        let zipped = exon_start_vec.iter().zip(exon_end_vec.iter());
        let mut tpos = 0; // track how many transcript bases has been seen
        for (exon_start, exon_end) in zipped {
            let exon = Exon::new(gstart, gend, tpos + 1);
            tpos = tpos + exon.exon_size;
            exon_vec.push(exon);
            }
        }
        exon_vec
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
        for exon in self.exons.iter() {
            if max(exon.gstart, &start) <= min(exon.gend, &end) {
                return Ok(true);
            }
        }
        return Ok(false);
    }

    /// ---
    /// 
    /// Get block sizes for in the genomic ranges for the transcript range
    /// 
    /// Args:
    ///     tstart (int): leftmost position on the transcript
    ///     tend (int): rightmost position on the transcript
    /// 
    /// return:
    ///     tuple(list, list): A tuple of list of block starts, and block ends
    pub fn blocks(&self, tstart: i32, tend: i32) -> PyResult<(Vec<i32>, Vec<i32>)>{
        assert!(tend <= self.transcript_length);
        assert!(tstart < tend);
        assert!(tstart > 0);
        let mut block_starts: Vec<i32> = vec![];
        let mut block_ends: Vec<i32> = vec![];
        let block_start: i32;
        let block_end: i32;
        let exon_start = 0;
        let exon_end = 0;
        let mut collect = 0;
        let mut start_offset = 0
        let mut end_offset = 0
        
        for exon in self.exons.iter() {
            if exon.tcontains(tstart) {
                start_offset = tstart - exon.tstart;
                if exon.tcontains(tend) {
                    /* contain CDS and CDE 
                    CDS           CDE
                        |->         ||
                    |===================|-----------|============|
                        (this exon)
                    */
                    collect = 1;
                    end_offset = tend - exon.tstart;
                    if self.strand == "-" {
                        block_start = exon.gend - end_offset;
                        block_end = exon.gend - start_offset;
                    } else {
                        block_start = exon.gstart + start_offset;
                        block_end = exon.gstart + end_offset;
                    }
                }
                else {


                }

            } else if collect != 0{

            }


                    if exon_start <= cds_end && cds_end <= exon_end {
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





        }

    }
        Ok((block_starts, block_ends))
}

fn adding(starts: Vec<i32>, sizes: Vec<i32>) -> (Vec<i32>,Vec<i32>){
    let mut res = vec![];
    for (start, size) in starts.iter().zip(sizes.iter()){
        let sum = start + size;
        res.push(sum);
    }
    (starts, res)
}