#![allow(dead_code)]

use std::collections::BTreeMap;
use std::collections::HashSet;

use crate::Feature;
use crate::Interval;
use crate::Node;

pub struct IntervalTree {
    pub all_intervals: HashSet<Interval>,
    pub top_node: Option<Box<Node>>,
    pub boundary_table: BTreeMap<i32,i32>
}

impl IntervalTree {
    pub fn new(intervals: Option<Vec<Interval>>) -> Self {
        let Some(mut intervals) = intervals else {
            return IntervalTree {
                all_intervals: HashSet::new(),
                top_node: None,
                boundary_table: BTreeMap::new(),
            };
        };

        if intervals.is_empty() {
            return IntervalTree {
                all_intervals: HashSet::new(),
                top_node: None,
                boundary_table: BTreeMap::new(),
            };
        }

        // Group equal feature IDs/strands next to each other. This replaces
        // the old HashSet -> HashMap<Vec<Interval>> -> HashSet pipeline and
        // lets us merge in one linear pass after sorting.
        intervals.sort_unstable_by(|a, b| {
            let af = a.data.as_ref().unwrap();
            let bf = b.data.as_ref().unwrap();
            af.name()
                .cmp(bf.name())
                .then(af.strand().cmp(&bf.strand()))
                .then(a.start.cmp(&b.start))
                .then(a.end.cmp(&b.end))
        });
        intervals.dedup();

        let mut merged: Vec<Interval> = Vec::with_capacity(intervals.len());
        for interval in intervals {
            let should_merge = merged.last().map(|last| {
                let last_feature = last.data.as_ref().unwrap();
                let feature = interval.data.as_ref().unwrap();
                last_feature.name() == feature.name()
                    && last_feature.strand() == feature.strand()
                    && last.end >= interval.start.saturating_sub(1)
            }).unwrap_or(false);

            if should_merge {
                let last = merged.last_mut().unwrap();
                if interval.end > last.end {
                    last.end = interval.end;
                }
            } else {
                merged.push(interval);
            }
        }

        // The static tree is queried by genomic position.
        merged.sort_unstable();

        let all_intervals: HashSet<Interval> = merged.iter().cloned().collect();
        let top_node = Node::from_sorted_slice(&merged);
        let mut it = IntervalTree {
            all_intervals,
            top_node,
            boundary_table: BTreeMap::new(),
        };

        // Keep the legacy boundary table for the older mutation/debug APIs,
        // although the counting overlap path no longer depends on it.
        for interval in &merged {
            it.add_boundaries(interval);
        }

        if let Some(root) = it.top_node.as_mut() {
            root.update_max_ends();
        }
        it
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

    pub fn overlap(&self, start: i32, end: i32) -> Vec<&Interval> {
        if self.is_empty() || start > end {
            return Vec::new();
        }

        let root = self.top_node.as_ref().unwrap();
        let mut result = Vec::new();
        root.search_overlap_range_into(start, end, &mut result);
        result
    }

    pub fn contains(&self, start: i32, end: i32) -> Vec<&Interval> {
        let mut result: Vec<&Interval> = self.overlap(start, end);
        result.retain(|x| x.start <= start && x.end >= end);
        result
    }

}