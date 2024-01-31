use crate::Feature;
use crate::Interval;

pub struct IntervalTree {
    intervals: Vec<Interval>,
}

impl IntervalTree {
    pub fn new() -> Self {
        IntervalTree {
            intervals: Vec::new(),
        }
    }

    pub fn insert(&mut self, start: i32, end: i32, data: Option<Feature>) {
        self.intervals.push(Interval::new(start, end, data));
    }
}