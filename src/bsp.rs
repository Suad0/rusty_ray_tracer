use std::vec::Vec;

#[derive(Debug, Clone)]
pub(crate) struct Point {
    pub(crate) x: f32,
    pub(crate) y: f32,
}

#[derive(Debug)]
pub struct BSPNode {
    points: Vec<Point>, // Field to store points within the node
    left: Option<Box<BSPNode>>,
    right: Option<Box<BSPNode>>,
    line_x: f32,
}

impl BSPNode {
    fn new(pts: Vec<Point>, line_x: f32) -> BSPNode {
        BSPNode {
            points: pts,  // Pass the provided points to the field
            left: None,
            right: None,
            line_x,
        }
    }

    fn add_point(&mut self, point: Point) { // Now takes point by reference
        if point.x < self.line_x {
            if self.left.is_none() {
                self.left = Some(Box::new(BSPNode::new(Vec::new(), self.line_x - 0.5)));
            }
            self.left.as_mut().unwrap().add_point(point); // Pass reference
        } else {
            if self.right.is_none() {
                self.right = Some(Box::new(BSPNode::new(Vec::new(), self.line_x + 0.5)));
            }
            self.right.as_mut().unwrap().add_point(point); // Pass reference
        }
    }

    fn print(&self) {
        if let Some(ref left) = self.left {
            left.print();
        }
        for point in &self.points { // Use self.points for iteration
            println!("({}, {})", point.x, point.y);
        }
        if let Some(ref right) = self.right {
            right.print();
        }
    }
}

#[derive(Debug)]
pub(crate) struct BSPTree {
    root: Option<Box<BSPNode>>,
}

impl BSPTree {
    pub(crate) fn new() -> BSPTree {
        BSPTree { root: None }
    }

    pub(crate) fn add_point(&mut self, point: Point) { // Now takes point by reference
        if self.root.is_none() {
            self.root = Some(Box::new(BSPNode::new(vec![point.clone()], point.x))); // Clone the point
        } else {
            self.root.as_mut().unwrap().add_point(point); // Pass reference
        }
    }

    pub(crate) fn print(&self) {
        if let Some(ref root) = self.root {
            root.print();
        }
        println!();
    }
}

fn main() {
    let mut tree = BSPTree::new();
    tree.add_point(Point { x: 5.0, y: 3.0 });
    tree.add_point(Point { x: 1.0, y: 7.0 });
    tree.add_point(Point { x: 10.0, y: 2.0 });
    tree.add_point(Point { x: 3.0, y: 6.0 });
    tree.add_point(Point { x: 8.0, y: 1.0 });

    println!("BSP Tree points:");
    tree.print();
}