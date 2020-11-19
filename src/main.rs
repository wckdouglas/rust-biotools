mod transcriptome;
mod utils;
use std::string::String;

fn main(){
    let txome = transcriptome::Transcriptome::new(String::from("test/data/test.refFlat"));
    println!("{}", txome.get("MIR6859-1").chrom)
}
