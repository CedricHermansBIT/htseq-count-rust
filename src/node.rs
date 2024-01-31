#![allow(dead_code)]

use crate::interval::Interval;
use std::collections::HashSet;

#[derive(Debug, Eq, PartialEq, Clone)]
pub struct Node {
    pub x_center: i32,
    pub s_center: HashSet<Interval>,
    left_node: Option<Box<Node>>,
    right_node: Option<Box<Node>>,
    pub depth: i32,
    pub balance: i32,
}

impl Node {
    pub fn new(x_center: i32, s_center: HashSet<Interval>) -> Self {
        Node {
            x_center,
            s_center,
            left_node: None,
            right_node: None,
            depth: 0,
            balance: 0,
        }
    }

    pub fn new_empty() -> Self {
        Node {
            x_center: 0,
            s_center: HashSet::new(),
            left_node: None,
            right_node: None,
            depth: 0,
            balance: 0,
        }
    }

    pub fn from_interval(interval: Interval) -> Self {
        let mut s_center = HashSet::new();
        s_center.insert(interval.clone());
        Node::new(interval.start, s_center)
    }

    pub fn from_intervals(intervals: Vec<Interval>) -> Option<Box<Node>> {
        if intervals.is_empty() {
            None
        } else {
            let mut sorted_intervals = intervals;
            sorted_intervals.sort();
            Node::from_sorted_intervals(sorted_intervals)
        }
    }

    pub fn from_sorted_intervals(intervals: Vec<Interval>) -> Option<Box<Node>> {
        if intervals.is_empty() {
            None
        } else {
            let mid = intervals.len() / 2;
            let mut s_center = HashSet::new();
            s_center.insert(intervals[mid].clone());

            let mut node = Node::new(intervals[mid].start, s_center);

            let left_intervals = intervals[0..mid].to_vec();
            let right_intervals = intervals[mid+1..].to_vec();

            node.left_node = Node::from_sorted_intervals(left_intervals);
            node.right_node = Node::from_sorted_intervals(right_intervals);

            let mut node = Some(Box::new(node));
            // rotate the tree to balance it
            node = node?.rotate();
            node

        }

    }

    pub fn center_hit(&self, interval: &Interval) -> bool {
        interval.contains_point(self.x_center)
    }

    pub fn hit_branch(&self, interval: &Interval) -> bool {
        interval.start > self.x_center
    }

    pub fn refresh_balance(&mut self) {
        let left_depth = self.left_node.as_ref().map_or(0, |n| n.depth);
        let right_depth = self.right_node.as_ref().map_or(0, |n| n.depth);
        self.depth = 1 + std::cmp::max(left_depth, right_depth);
        self.balance = right_depth - left_depth;
    }

    pub fn compute_depth(&mut self) -> i32 {
        let left_depth = self.left_node.as_mut().map_or(0, |n| n.compute_depth());
        let right_depth = self.right_node.as_mut().map_or(0, |n| n.compute_depth());
        self.depth = 1 + std::cmp::max(left_depth, right_depth);
        self.depth
    }

    pub fn rotate(&mut self) -> Option<Box<Node>> {
        self.refresh_balance();
        if self.balance.abs() < 2 {
            return Some(Box::new(self.clone()));
        }
        let my_heavy = self.balance > 0;
        let child_heavy = match my_heavy {
            true => self.right_node.as_ref()?.balance > 0,
            false => self.left_node.as_ref()?.balance > 0,
        };
        if my_heavy == child_heavy || self.right_node.as_ref()?.balance == 0 {
            self.srotate()
        } else {
            self.drotate()
        }
    }

    pub fn srotate(&mut self) -> Option<Box<Node>> {
        let heavy = self.balance > 0;
        let mut save = match heavy {
            true => self.right_node.take(),
            false => self.left_node.take(),
        }?;
        match heavy {
            true => self.right_node = save.left_node.take(),
            false => self.left_node = save.right_node.take(),
        };
        save.rotate();  // Needed to ensure the 2 and 3 are balanced under new subnode
        Some(save)
    }

    pub fn drotate(&mut self) -> Option<Box<Node>> {
        let my_heavy = self.balance > 0;
        match my_heavy {
            true => self.right_node.as_mut()?.srotate(),
            false => self.left_node.as_mut()?.srotate(),
        };
        self.refresh_balance();
        self.srotate()
    }


    pub fn add(&mut self, interval: Interval) -> &mut Self {
        if self.center_hit(&interval) {
            self.s_center.insert(interval);
        } else {
            let direction = self.hit_branch(&interval);
            match direction {
                true => {
                    if let Some(right_node) = &mut self.right_node {
                        right_node.add(interval);
                    } else {
                        self.right_node = Some(Box::new(Node::from_interval(interval)));
                    }
                },
                false => {
                    if let Some(left_node) = &mut self.left_node {
                        left_node.add(interval);
                    } else {
                        self.left_node = Some(Box::new(Node::from_interval(interval)));
                    }
                },
            }
            self.refresh_balance();
        }
        self
    }

    pub fn remove(&mut self, interval: Interval) -> Result<&mut Self, &'static str> {
        self.remove_interval_helper(interval, true)
    }

    pub fn discard(&mut self, interval: Interval) -> Result<&mut Self, &'static str> {
        self.remove_interval_helper(interval, false)
    }

    pub fn remove_interval_helper(&mut self, interval: Interval, should_raise_error: bool) -> Result<&mut Self, &'static str> {
        if self.center_hit(&interval) {
            if !self.s_center.remove(&interval) && should_raise_error {
                return Err("Interval not found");
            }
            if self.s_center.is_empty() {
                // implement prune method here
            }
        } else {
            let direction = self.hit_branch(&interval);
            match direction {
                true => {
                    if let Some(right_node) = &mut self.right_node {
                        right_node.remove_interval_helper(interval, should_raise_error)?;
                    } else if should_raise_error {
                        return Err("Interval not found");
                    }
                },
                false => {
                    if let Some(left_node) = &mut self.left_node {
                        left_node.remove_interval_helper(interval, should_raise_error)?;
                    } else if should_raise_error {
                        return Err("Interval not found");
                    }
                },
            }
            self.refresh_balance();
        }
        Ok(self)
    }

    pub fn search_overlap(&self, point_list: Vec<i32>) -> Vec<Interval> {
        let mut result = Vec::new();
        for point in point_list {
            result.extend(self.search_point(point));
        }
        result
    }

    pub fn search_point(&self, point: i32) -> Vec<Interval> {
        let mut result = Vec::new();
        for interval in &self.s_center {
            if interval.contains_point(point) {
                result.push(interval.clone());
            }
        }
        if point < self.x_center {
            if let Some(left_node) = &self.left_node {
                result.extend(left_node.search_point(point));
            }
        } else if point > self.x_center {
            if let Some(right_node) = &self.right_node {
                result.extend(right_node.search_point(point));
            }
        }
        result
    }

    pub fn contains_point(&self, point: i32) -> bool {
        for interval in &self.s_center {
            if interval.contains_point(point) {
                return true;
            }
        }
        if point < self.x_center {
            if let Some(left_node) = &self.left_node {
                return left_node.contains_point(point);
            }
        } else if point > self.x_center {
            if let Some(right_node) = &self.right_node {
                return right_node.contains_point(point);
            }
        }
        false
    }

    pub fn all_children(&self) -> Vec<Interval> {
        let mut result = Vec::new();
        for interval in &self.s_center {
            result.push(interval.clone());
        }
        if let Some(left_node) = &self.left_node {
            result.extend(left_node.all_children());
        }
        if let Some(right_node) = &self.right_node {
            result.extend(right_node.all_children());
        }
        result
    }

    pub fn count_nodes(&self) -> usize {
        let mut result = 1;
        if let Some(left_node) = &self.left_node {
            result += left_node.count_nodes();
        }
        if let Some(right_node) = &self.right_node {
            result += right_node.count_nodes();
        }
        result
    }
    // indent by default is 0
    pub fn print_structure(&self, indent: usize, to_string:bool) -> Vec<String> {
        let newline = "\n";
        let spaces = " ".repeat(indent);

        let mut result = vec![format!("{}{}", self, newline)];
        if self.s_center.len() > 0 {
            for interval in &self.s_center {
                result.push(format!("{} {}{}", spaces, interval, newline));
            }
        }
        if let Some(left_node) = &self.left_node {
            result.push(format!("{}{}{}", spaces, "Left:", newline));
            result.extend(left_node.print_structure(indent+1,true));
        }
        if let Some(right_node) = &self.right_node {
            result.push(format!("{}{}{}", spaces, "Right:", newline));
            result.extend(right_node.print_structure(indent+1,true));
        }
        
        // join the vector into a string
        let result = result.join("");

        if to_string {
            return vec![result]
        } else {
            println!("{}", result);
            return vec![]
        }

    }
}
    
// display format

impl std::fmt::Display for Node {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        //"Node<{0}, depth={1}, balance={2}>"
        write!(f, "Node<{}, depth={}, balance={}>", self.x_center, self.depth, self.balance)
    }
}