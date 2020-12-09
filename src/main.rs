mod transcriptome;
mod utils;
use std::string::String;

fn main(){
    let mut txome = transcriptome::Transcriptome(String::from("test/data/test.refFlat"));
    let gene_name = String::from("SNX3");
    //let gene = txome.fetch_gene(&gene_name);
    //let transcript = &gene.transcripts[0].rnaId;
    //println!("{}", transcript);
}
