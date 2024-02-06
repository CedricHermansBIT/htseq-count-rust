// Struct to store the features
#[derive(Debug, Eq, PartialEq, Hash, Clone, Default)]
pub struct Feature {
    //type_: String,
    pub name: String,
    pub chr: i32,
    pub start: i32,
    pub end: i32,
}