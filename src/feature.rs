// Struct to store the features
#[derive(Debug, Eq, PartialEq, Hash, Clone)]
pub struct Feature {
    //type_: String,
    pub name: String,
    pub chr: i32,
    pub start: i32,
    pub end: i32,
}

// default implementation for Feature
impl Default for Feature {
    fn default() -> Self {
        Feature {
            name: String::from(""),
            chr: 0,
            start: 0,
            end: 0
        }
    }
}