// Struct to store the features
#[derive(Debug, Eq, PartialEq, Hash, Clone)]
pub struct Feature {
    //type_: String,
    pub name: String,
    pub chr: String,
    pub start: u32,
    pub end: u32,
}

// default implementation for Feature
impl Default for Feature {
    fn default() -> Self {
        Feature {
            name: String::from(""),
            chr: String::from(""),
            start: 0,
            end: 0
        }
    }
}