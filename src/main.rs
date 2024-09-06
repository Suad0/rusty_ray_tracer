use std::ops::{Add, Div, Mul, Sub};
use rand::Rng;

#[derive(Copy, Clone, Debug)]
struct Vec3 {
    x: f64,
    y: f64,
    z: f64,
}

impl Mul<f64> for Vec3 {
    type Output = Vec3;

    fn mul(self, rhs: f64) -> Vec3 {
        Vec3 {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}


impl Div<f64> for &Vec3 {
    type Output = Vec3;

    fn div(self, rhs: f64) -> Self::Output {
        Vec3 {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}

impl Div<f64> for Vec3 {
    type Output = Vec3;

    fn div(self, rhs: f64) -> Self::Output {
        Vec3 {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}


impl Add<Vec3> for Vec3 {
    type Output = Vec3;

    fn add(self, rhs: Vec3) -> Vec3 {
        Vec3 {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl Sub for Vec3 {
    type Output = Vec3;

    fn sub(self, rhs: Self) -> Self::Output {
        Vec3 {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl Vec3 {
    fn new(x: f64, y: f64, z: f64) -> Vec3 {
        Vec3 { x, y, z }
    }

    fn dot(self, rhs: Self) -> f64 {
        self.x * rhs.x + self.y * rhs.y + self.z * rhs.z
    }

    fn unit_vector(&self) -> Vec3 {
        let length = self.length();
        if length == 0.0 {
            panic!("Cannot normalize a zero vector");
        }
        self / length
    }

    fn length(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }
}

struct Ray {
    origin: Vec3,
    direction: Vec3,
}

impl Ray {
    fn new(origin: Vec3, direction: Vec3) -> Ray {
        Ray { origin, direction }
    }

    fn point_at_parameter(&self, t: f64) -> Vec3 {
        self.origin + self.direction * t
    }
}

struct Sphere {
    center: Vec3,
    radius: f64,
}

impl Sphere {
    fn new(center: Vec3, radius: f64) -> Sphere {
        Sphere { center, radius }
    }

    fn hit(&self, ray: &Ray) -> Option<f64> {
        let oc = ray.origin - self.center;
        let a = ray.direction.dot(ray.direction);
        let b = 2.0 * oc.dot(ray.direction);
        let c = oc.dot(oc) - self.radius * self.radius;
        let discriminant = b * b - 4.0 * a * c;
        if discriminant > 0.0 {
            Some((-b - discriminant.sqrt()) / (2.0 * a))
        } else {
            None
        }
    }
}

struct Plane {
    point: Vec3,    // A point on the plane
    normal: Vec3,   // The plane's normal vector
}

impl Plane {
    fn new(point: Vec3, normal: Vec3) -> Plane {
        Plane { point, normal: normal.unit_vector() }
    }

    fn hit(&self, ray: &Ray) -> Option<f64> {
        let denom = self.normal.dot(ray.direction);
        if denom.abs() > 1e-6 {
            let t = (self.point - ray.origin).dot(self.normal) / denom;
            if t >= 0.0 {
                return Some(t);
            }
        }
        None
    }
}

struct Camera {
    origin: Vec3,
    lower_left_corner: Vec3,
    horizontal: Vec3,
    vertical: Vec3,
}

impl Camera {
    fn new() -> Camera {
        Camera {
            origin: Vec3::new(0.0, 0.0, 0.0),
            lower_left_corner: Vec3::new(-2.0, -1.0, -1.0),
            horizontal: Vec3::new(4.0, 0.0, 0.0),
            vertical: Vec3::new(0.0, 2.0, 0.0),
        }
    }

    fn get_ray(&self, u: f64, v: f64) -> Ray {
        Ray::new(
            self.origin,
            self.lower_left_corner + self.horizontal * u + self.vertical * v - self.origin,
        )
    }
}

struct Light {
    position: Vec3,
    intensity: f64, // Brightness of the light
}

impl Light {
    fn new(position: Vec3, intensity: f64) -> Light {
        Light { position, intensity }
    }
}

fn in_shadow(ray: &Ray, spheres: &[Sphere]) -> bool {
    for sphere in spheres {
        if let Some(_) = sphere.hit(ray) {
            return true; // There is an object blocking the light
        }
    }
    false
}

fn color(ray: &Ray, spheres: &[Sphere], planes: &[Plane], light: &Light) -> Vec3 {
    let mut closest_t = f64::MAX;
    let mut hit_color = None;

    // Check sphere intersections
    for sphere in spheres {
        if let Some(t) = sphere.hit(ray) {
            if t < closest_t {
                closest_t = t;
                let hit_point = ray.point_at_parameter(t);
                let N = (hit_point - sphere.center).unit_vector();

                // Check for shadows
                let light_dir = (light.position - hit_point).unit_vector();
                let shadow_ray = Ray::new(hit_point, light_dir);
                if in_shadow(&shadow_ray, spheres) {
                    hit_color = Some(Vec3::new(0.1, 0.1, 0.1)); // In shadow, dark color
                } else {
                    // Diffuse lighting (Lambertian reflection)
                    let light_intensity = light.intensity * N.dot(light_dir).max(0.0);
                    hit_color = Some(Vec3::new(1.0, 0.5, 0.5) * light_intensity); // Diffuse shading
                }
            }
        }
    }

    // Check plane intersections
    for plane in planes {
        if let Some(t) = plane.hit(ray) {
            if t < closest_t {
                closest_t = t;
                let hit_point = ray.point_at_parameter(t);
                let N = plane.normal;

                // Check for shadows
                let light_dir = (light.position - hit_point).unit_vector();
                let shadow_ray = Ray::new(hit_point, light_dir);
                if in_shadow(&shadow_ray, spheres) {
                    hit_color = Some(Vec3::new(0.1, 0.1, 0.1)); // In shadow, dark color
                } else {
                    // Diffuse lighting (Lambertian reflection)
                    let light_intensity = light.intensity * N.dot(light_dir).max(0.0);
                    hit_color = Some(Vec3::new(0.5, 0.5, 1.0) * light_intensity); // Diffuse shading
                }
            }
        }
    }

    // Return hit color, or a default sky color if no object was hit
    if let Some(color) = hit_color {
        color
    } else {
        let unit_direction = ray.direction.unit_vector();
        let t = 0.5 * (unit_direction.y + 1.0);
        Vec3::new(1.0 - t, 1.0 - t, 1.0 - t) + Vec3::new(t * 0.5, t * 0.7, t)
    }
}

fn main() {
    let nx = 200;
    let ny = 100;
    let ns = 10; // Number of samples per pixel for anti-aliasing

    println!("P3\n{} {}\n255", nx, ny);
    let camera = Camera::new();
    let light = Light::new(Vec3::new(5.0, 5.0, 2.0), 1.0); // Point light

    let spheres = vec![
        Sphere::new(Vec3::new(0.0, 0.0, -1.0), 0.5),
        Sphere::new(Vec3::new(1.0, 0.0, -1.0), 0.5),
        Sphere::new(Vec3::new(-1.0, 0.0, -1.0), 0.5),
    ];

    let planes = vec![
        Plane::new(Vec3::new(0.0, -0.5, 0.0), Vec3::new(0.0, 1.0, 0.0)), // Ground plane
    ];

    let mut rng = rand::thread_rng();

    for j in (0..ny).rev() {
        for i in 0..nx {
            let mut col = Vec3::new(0.0, 0.0, 0.0);
            for _ in 0..ns {
                let u = (i as f64 + rng.gen::<f64>()) / nx as f64;
                let v = (j as f64 + rng.gen::<f64>()) / ny as f64;
                let ray = camera.get_ray(u, v);
                col = col + color(&ray, &spheres, &planes, &light);
            }
            col = col / ns as f64;
            let ir = (255.99 * col.x) as i32;
            let ig = (255.99 * col.y) as i32;
            let ib = (255.99 * col.z) as i32;
            println!("{} {} {}", ir, ig, ib);
        }
    }
}
