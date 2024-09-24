use crate::ray_tracer::{Plane, Ray, Sphere};
use std::any::Any;


#[derive(Debug, Clone)]
struct BSPNode {
    objects: Vec<Box<dyn Hittable>>,  // Objects (like spheres, planes)
    left: Option<Box<BSPNode>>,
    right: Option<Box<BSPNode>>,
    split_axis: char,  // 'x', 'y', or 'z'
    split_position: f64,  // Split position along the axis
}

// Trait for objects that can be hit by a ray
trait Hittable {
    fn hit(&self, ray: &Ray) -> Option<f64>;
    fn as_any(&self) -> &dyn Any;
}

impl Hittable for Sphere {
    fn hit(&self, ray: &Ray) -> Option<f64> {
        self.hit(ray)
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}

impl Hittable for Plane {
    fn hit(&self, ray: &Ray) -> Option<f64> {
        self.hit(ray)
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}


fn get_object_position(object: &Box<dyn Hittable>, axis: char) -> f64 {
    // Downcast the Hittable trait object to the concrete types (Sphere, Plane, etc.)
    if let Some(sphere) = object.as_any().downcast_ref::<Sphere>() {
        match axis {
            'x' => sphere.center.x,
            'y' => sphere.center.y,
            'z' => sphere.center.z,
            _ => panic!("Invalid axis"),
        }
    } else if let Some(plane) = object.as_any().downcast_ref::<Plane>() {
        match axis {
            'x' => plane.point.x,
            'y' => plane.point.y,
            'z' => plane.point.z,
            _ => panic!("Invalid axis"),
        }
    } else {
        panic!("Unsupported object type in BSP");
    }
}

fn build_bsp(objects: Vec<Box<dyn Hittable>>, depth: u32) -> BSPNode {
    if objects.len() == 1 {
        return BSPNode {
            objects,
            left: None,
            right: None,
            split_axis: 'x',
            split_position: 0.0,
        };
    }

    // Choose axis based on depth
    let axis = match depth % 3 {
        0 => 'x',
        1 => 'y',
        _ => 'z',
    };

    // Sort objects along the chosen axis
    let mut sorted_objects = objects.clone();
    sorted_objects.sort_by(|a, b| {
        let a_pos = get_object_position(a, axis);
        let b_pos = get_object_position(b, axis);
        a_pos.partial_cmp(&b_pos).unwrap()
    });

    // Split the objects at the median
    let mid = sorted_objects.len() / 2;
    let left_objects = sorted_objects[..mid].to_vec();
    let right_objects = sorted_objects[mid..].to_vec();

    BSPNode {
        objects: vec![],
        left: Some(Box::new(build_bsp(left_objects, depth + 1))),
        right: Some(Box::new(build_bsp(right_objects, depth + 1))),
        split_axis: axis,
        split_position: get_object_position(&sorted_objects[mid], axis),
    }
}








