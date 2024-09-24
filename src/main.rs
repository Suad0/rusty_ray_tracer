use crate::bsp::{BSPTree, Point};

mod bsp;
mod ray_tracer;
mod bsp_node;


fn main() {
    let mut tree = BSPTree::new(5);  // Set max depth to 5 to prevent infinite recursion
    tree.add_point(Point { x: 5.0, y: 3.0 });
    tree.add_point(Point { x: 1.0, y: 7.0 });
    tree.add_point(Point { x: 10.0, y: 2.0 });
    tree.add_point(Point { x: 3.0, y: 6.0 });
    tree.add_point(Point { x: 8.0, y: 1.0 });

    println!("BSP Tree points:");
    tree.print();
}
