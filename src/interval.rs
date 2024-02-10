#![allow(dead_code)]

use crate::Feature;

#[derive(Debug, Eq, PartialEq, Hash, Clone)]
pub struct Interval {
    pub start: i32,
    pub end: i32,
    //optional data
    pub data: Option<Feature>,
}

impl Interval {
    pub fn new(start: i32, end: i32, data: Option<Feature>) -> Self {
        Interval { start, end, data }
    }

    pub fn overlaps(&self, other: &Interval) -> bool {
        other.start < self.end && other.end > self.start
    }

    pub fn overlap_size(&self, other: &Interval) -> i32 {
        if self.overlaps(other) {
            let start = std::cmp::max(self.start, other.start);
            let end = std::cmp::min(self.end, other.end);
            end - start
        } else {
            0
        }
    }

    pub fn contains_point(&self, point: i32) -> bool {
        self.start <= point && point <= self.end
    }

    pub fn range_matches(&self, other: &Interval) -> bool {
        self.start == other.start && self.end == other.end
    }

    pub fn contains_interval(&self, other: &Interval) -> bool {
        self.start <= other.start && self.end >= other.end
    }

    pub fn distance_to(&self, other: &Interval) -> i32 {
        if self.overlaps(other) {
            0
        } else if self.start > other.end {
            self.start - other.end
        } else {
            other.start - self.end
        }
        
    }

    pub fn distance_to_point(&self, point: i32) -> i32 {
        if self.contains_point(point) {
            0
        } else if self.start > point {
            self.start - point
        } else {
            point - self.end
        }
    }

    pub fn is_null(&self) -> bool {
        self.start>=self.end
    }

    pub fn length(&self) -> i32 {
        if self.is_null() {
            0
        } else {
            self.end - self.start
        }
    }

    fn raise_if_null(&self, other: &Interval) {
        if self.is_null() || other.is_null(){
            panic!("Cannot perform operation on null interval");
        }
    }
    pub fn lt(&self, other: &Interval) -> bool {
        // Strictly less than
        self.raise_if_null(other);
        self.end <= other.start
    }

    pub fn le(&self, other: &Interval) -> bool {
        // Less than or overlaps
        self.raise_if_null(other);
        self.end < other.end
    }

    pub fn gt(&self, other: &Interval) -> bool {
        // Strictly greater than
        self.raise_if_null(other);
        self.start > other.end
    }

    pub fn ge(&self, other: &Interval) -> bool {
        // Greater than or overlaps
        self.raise_if_null(other);
        self.start >= other.start
    }

    pub fn overlaps_with(&self, start: i32, end: i32) -> bool {
        self.start < end && self.end-1 > start
    }

}

impl PartialOrd for Interval {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Interval {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // check if both intervals are not null
        // if self.is_null() || other.is_null() {
        //     // throw error
        //     panic!("Cannot compare null intervals");
        // }

        // first by start, then by end, then by data.name (alphabetically)
        self.start.cmp(&other.start)
            .then(self.end.cmp(&other.end))
            .then(self.data.as_ref().unwrap().name().cmp(&other.data.as_ref().unwrap().name()))

    }
}

impl std::fmt::Display for Interval {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match &self.data {
            Some(feature) => write!(f, "Interval({}, {}, {})", self.start, self.end, feature.name()),
            None => write!(f, "Interval({}, {})", self.start, self.end),
        }
    }
}


