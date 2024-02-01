#![allow(dead_code)]

use std::collections::HashMap;
use std::collections::HashSet;

use crate::Feature;
use crate::Interval;
use crate::Node;

pub struct IntervalTree {
    pub all_intervals: HashSet<Interval>,
    pub top_node: Option<Box<Node>>,
    pub boundary_table: HashMap<i32,i32>
}

impl IntervalTree {
    pub fn new(intervals: Option<Vec<Interval>>) -> Self {
        if let Some(intervals) = intervals {
            eprintln!("creating IntervalTree");
            // get unique intervals
            let mut unique_intervals: HashSet<Interval> = HashSet::new();
            for interval in intervals {
                unique_intervals.insert(interval);
            }
            // convert to vector
            //let unique_intervals: Vec<Interval> = unique_intervals.into_iter().collect();
            eprintln!("filtered intervals");
            // check if intervals are not null
            for interval in &unique_intervals {
                if interval.is_null() {
                    panic!("IntervalTree: Null Interval objects are not allowed in IntervalTree");
                }
            }
            let mut it = IntervalTree {
                all_intervals: unique_intervals.clone(),
                top_node: Node::from_intervals(unique_intervals.clone()),
                boundary_table: HashMap::new(),
            };
            eprintln!("created top node");

            for interval in &unique_intervals {
                it.add_boundaries(interval);
            }
            eprintln!("added boundaries");
            return it;
        }
        else {
            IntervalTree {
                all_intervals: HashSet::new(),
                top_node: None,
                boundary_table: HashMap::new(),
            }
        }
    }

    pub fn from_tuples(tuples: Vec<(i32,i32, Option<Feature>)>) -> IntervalTree {
        let mut intervals: Vec<Interval> = Vec::new();
        eprintln!("creating intervals");
        for tuple in tuples {
            let interval = Interval::new(tuple.0, tuple.1, tuple.2);
            intervals.push(interval);
        }
        eprintln!("created intervals");
        IntervalTree::new(Some(intervals))
    }

    fn add_boundaries(&mut self, interval: &Interval) {
        let start = interval.start;
        let end = interval.end;
        *self.boundary_table.entry(start).or_insert(0) += 1;
        *self.boundary_table.entry(end).or_insert(0) += 1;
    }

    fn remove_boundaries(&mut self, interval: &Interval) {
        let start = interval.start;
        let end = interval.end;
        if let Some(count) = self.boundary_table.get_mut(&start) {
            *count -= 1;
            if *count == 0 {
                self.boundary_table.remove(&start);
            }
        }
        if let Some(count) = self.boundary_table.get_mut(&end) {
            *count -= 1;
            if *count == 0 {
                self.boundary_table.remove(&end);
            }
        }
    }

    pub fn add(&mut self, interval: Interval) {
        if interval.is_null() {
            panic!("IntervalTree: Null Interval objects are not allowed in IntervalTree");
        }
        // check if interval is unique
        if self.all_intervals.contains(&interval) {
            return;
        }
        
        self.all_intervals.insert(interval.clone());
        self.add_boundaries(&interval);
        if let Some(top_node) = &mut self.top_node {
            top_node.add(interval);
        } else {
            let mut intervals: HashSet<Interval> = HashSet::new();
            intervals.insert(interval.clone());
            self.top_node = Node::from_intervals(intervals);
        }
    }

    pub fn addi(&mut self, start: i32, end: i32, data: Option<Feature>) {
        let interval = Interval::new(start, end, data);
        self.add(interval);
    }

    pub fn update(&mut self, intervals: HashSet<Interval>) {
        for interval in intervals {
            self.add(interval);
        }
    }

    pub fn remove(&mut self, interval: &Interval) {
        if !self.all_intervals.contains(interval) {
            return;
        }
        self.all_intervals.retain(|x| !x.range_matches(interval));
        self.remove_boundaries(interval);
        if let Some(top_node) = &mut self.top_node {
            let _ = top_node.remove(interval.clone());
        }
    }

    pub fn removei(&mut self, start: i32, end: i32, data: Option<Feature>) {
        let interval = Interval::new(start, end, data);
        self.remove(&interval);
    }

    pub fn discard(&mut self, interval: &Interval) {
        if !self.all_intervals.contains(interval) {
            return;
        }
        self.all_intervals.retain(|x| !x.range_matches(interval));
        self.remove_boundaries(interval);
        if let Some(top_node) = &mut self.top_node {
            let _ = top_node.discard(interval.clone());
        }
    }

    pub fn discardi(&mut self, start: i32, end: i32, data: Option<Feature>) {
        let interval = Interval::new(start, end, data);
        self.discard(&interval);
    }

    pub fn difference(&self, other: &IntervalTree) -> IntervalTree {
        let mut result: Vec<Interval> = Vec::new();
        for interval in &self.all_intervals {
            if !other.all_intervals.contains(interval) {
                result.push(interval.clone());
            }
        }
        IntervalTree::new(Some(result))
    }

    pub fn difference_update(&mut self, other: &IntervalTree) {
        // discard intervals in self that are also in other
        for interval in &other.all_intervals {
            if self.all_intervals.contains(interval) {
                self.discard(interval);
            }
        }
    }

    pub fn union(&self, other: &IntervalTree) -> IntervalTree {
        // return a new IntervalTree with all intervals from self and other
        let mut result: Vec<Interval> = Vec::new();
        for interval in &self.all_intervals {
            result.push(interval.clone());
        }
        for interval in &other.all_intervals {
            if !result.contains(interval) {
                result.push(interval.clone());
            }
        }
        IntervalTree::new(Some(result))
    }

    pub fn intersection(&self, other: &IntervalTree) -> IntervalTree {
        // return a new IntervalTree with all intervals that are in both self and other
        let mut result: Vec<Interval> = Vec::new();
        let (shorter, longer) = if self.len() < other.len() {
            (&self.all_intervals, &other.all_intervals)
        } else {
            (&other.all_intervals, &self.all_intervals)
        };
        for interval in shorter {
            if longer.contains(interval) {
                result.push(interval.clone());
            }
        }
        IntervalTree::new(Some(result))
    }

    pub fn intersection_update(&mut self, other: &IntervalTree) {
        // remove intervals from self unless they are also in other
        let intervals = self.all_intervals.clone();
        for interval in intervals {
            if !other.all_intervals.contains(&interval) {
                self.remove(&interval);
            }
        }
    }

    pub fn symmetric_difference(&self, other: &IntervalTree) -> IntervalTree {
        // return a new IntervalTree with all intervals that are in either self or other but not both
        let mut result: Vec<Interval> = Vec::new();
        for interval in &self.all_intervals {
            if !other.all_intervals.contains(interval) {
                result.push(interval.clone());
            }
        }
        for interval in &other.all_intervals {
            if !self.all_intervals.contains(interval) {
                result.push(interval.clone());
            }
        }
        IntervalTree::new(Some(result))
    }

    pub fn symmetric_difference_update(&mut self, other: &mut IntervalTree) {
        // remove intervals from self unless they are also in other
        let intervals = self.all_intervals.clone();
        let other_intervals = other.all_intervals.clone();
        for interval in intervals {
            if other_intervals.contains(&interval) {
                self.remove(&interval);
                other.remove(&interval);
            }
        }
        self.update(other.all_intervals.clone());
    }

    pub fn clear(&mut self) {
        self.all_intervals.clear();
        self.top_node = None;
        self.boundary_table.clear();
    }

    pub fn is_empty(&self) -> bool {
        self.all_intervals.is_empty()
    }

    pub fn len(&self) -> usize {
        self.all_intervals.len()
    }

    pub fn overlaps_point(&self, point: i32) -> bool {
        if self.is_empty() {
            return false;
        }
        // use top_node to find overlapping intervals
        if let Some(top_node) = &self.top_node {
            return top_node.as_ref().contains_point(point)
        }
        false
    }

    pub fn overlaps_range(&self, start: i32, end: i32) -> bool {
        if self.is_empty() {
            return false;
        }
        if start >= end {
            return false;
        }
        if self.overlaps_point(start) {
            return true;
        }
        // begin < boundary < end
        for boundary in self.boundary_table.keys() {
            if start < *boundary && *boundary < end {
                return true;
            }
        }
        false
    }


}