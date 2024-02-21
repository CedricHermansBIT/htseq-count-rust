#![allow(dead_code)]

// Struct to store the features
#[derive(Debug, Eq, PartialEq, Hash, Clone, Default)]
pub struct Feature {
    //type_: String,
    name: String,
    chr: i32,
    start: i32,
    end: i32,
    strand: char,
}

impl Feature {
    pub fn new(name: String, chr: i32, start: i32, end: i32, strand: char) -> Self {
        Feature {
            name,
            chr,
            start,
            end,
            strand,
        }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn chr(&self) -> i32 {
        self.chr
    }

    pub fn start(&self) -> i32 {
        self.start
    }

    pub fn end(&self) -> i32 {
        self.end
    }

    pub fn strand(&self) -> char {
        self.strand
    }

    pub fn set_end(&mut self, end: i32) {
        self.end = end;
    }

}