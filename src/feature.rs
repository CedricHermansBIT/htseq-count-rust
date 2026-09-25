#![allow(dead_code)]

// Struct to store the features
#[derive(Debug, Eq, PartialEq, Hash, Clone, Default)]
pub struct Feature {
    // Integer ID is interned once during annotation loading. It is used in
    // per-read ambiguity logic so the hot path does not compare gene strings.
    id: usize,
    name: String,
    chr: i32,
    start: i32,
    end: i32,
    strand: char,
}

impl Feature {
    pub fn new(id: usize, name: String, chr: i32, start: i32, end: i32, strand: char) -> Self {
        Feature {
            id,
            name,
            chr,
            start,
            end,
            strand,
        }
    }

    pub fn id(&self) -> usize {
        self.id
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