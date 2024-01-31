// Struct to store the features
#[derive(Debug, Eq, PartialEq, Hash, Clone)]
pub struct Feature {
    //type_: String,
    pub name: String,
    pub chr: String,
    pub start: u32,
    pub end: u32,
    pub start_sorted_index: usize,
    pub end_sorted_index: usize,

}