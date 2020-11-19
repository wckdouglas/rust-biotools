mod transcriptome;
mod utils;
use std::string::String;

fn main(){
    transcriptome::Transcriptome::new(String::from("test/data/test.refFlat"));
}
