use std::vec::Vec;

#[derive(Debug, Clone)]
pub struct Point {
    pub x: f32,
    pub y: f32,
}

#[derive(Debug)]
pub struct BSPNode {
    points: Vec<Point>,  // Field to store points within the node
    left: Option<Box<BSPNode>>,
    right: Option<Box<BSPNode>>,
    line_x: f32,
    max_depth: usize,   // Add a max depth limit
    depth: usize,       // Track the depth of the current node
}

impl BSPNode {
    fn new(pts: Vec<Point>, line_x: f32, depth: usize, max_depth: usize) -> BSPNode {
        BSPNode {
            points: pts,
            left: None,
            right: None,
            line_x,
            depth,
            max_depth,
        }
    }

    fn add_point(&mut self, point: Point) {
        if self.depth >= self.max_depth || (point.x - self.line_x).abs() < 0.01 {
            // If we've reached max depth or the point is very close to the line, store it in this node
            self.points.push(point);
        } else if point.x < self.line_x {
            // Point goes to the left child
            if self.left.is_none() {
                self.left = Some(Box::new(BSPNode::new(Vec::new(), self.line_x - 0.5, self.depth + 1, self.max_depth)));
            }
            self.left.as_mut().unwrap().add_point(point);
        } else {
            // Point goes to the right child
            if self.right.is_none() {
                self.right = Some(Box::new(BSPNode::new(Vec::new(), self.line_x + 0.5, self.depth + 1, self.max_depth)));
            }
            self.right.as_mut().unwrap().add_point(point);
        }
    }

    fn print(&self) {
        if let Some(ref left) = self.left {
            left.print();
        }
        for point in &self.points {
            println!("({}, {})", point.x, point.y);
        }
        if let Some(ref right) = self.right {
            right.print();
        }
    }
}

#[derive(Debug)]
pub struct BSPTree {
    root: Option<Box<BSPNode>>,
    max_depth: usize,  // Add a max depth parameter
}

impl BSPTree {
    pub fn new(max_depth: usize) -> BSPTree {
        BSPTree { root: None, max_depth }
    }

    pub fn add_point(&mut self, point: Point) {
        if self.root.is_none() {
            self.root = Some(Box::new(BSPNode::new(vec![point.clone()], point.x, 0, self.max_depth)));
        } else {
            self.root.as_mut().unwrap().add_point(point);
        }
    }

    pub fn print(&self) {
        if let Some(ref root) = self.root {
            root.print();
        }
        println!();
    }
}


