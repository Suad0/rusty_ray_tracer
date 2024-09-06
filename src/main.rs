use std::ops::{Add, Div, Mul, Sub};
use image::{Rgb, RgbImage};
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
    lens_radius: f64,  // For depth of field
}

fn random_in_unit_disk() -> Vec3 {
    let mut rng = rand::thread_rng();
    loop {
        let p = Vec3::new(rng.gen::<f64>(), rng.gen::<f64>(), 0.0) * 2.0 - Vec3::new(1.0, 1.0, 0.0);
        if p.dot(p) < 1.0 {
            return p;
        }
    }
}

impl Camera {
    fn new(fov: f64, aspect: f64, aperture: f64, focus_dist: f64) -> Camera {
        let lens_radius = aperture / 2.0;
        let theta = fov.to_radians();
        let half_height = (theta / 2.0).tan();
        let half_width = aspect * half_height;
        let origin = Vec3::new(0.0, 0.0, 0.0);

        Camera {
            origin,
            lower_left_corner: origin - Vec3::new(half_width, half_height, 1.0) * focus_dist,
            horizontal: Vec3::new(2.0 * half_width, 0.0, 0.0) * focus_dist,
            vertical: Vec3::new(0.0, 2.0 * half_height, 0.0) * focus_dist,
            lens_radius,
        }
    }

    fn get_ray(&self, u: f64, v: f64) -> Ray {
        let rd = random_in_unit_disk() * self.lens_radius;
        let offset = rd;  // Modify this to implement depth of field
        Ray::new(
            self.origin + offset,
            self.lower_left_corner + self.horizontal * u + self.vertical * v - self.origin - offset,
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

fn in_shadow(ray: &Ray, spheres: &[Sphere], light: &Light) -> bool {
    let shadow_samples = 10;  // Number of random light samples
    let mut rng = rand::thread_rng();

    for _ in 0..shadow_samples {
        // Randomly jitter the light source to simulate area light
        let light_jitter = Vec3::new(
            rng.gen::<f64>() - 0.5,
            rng.gen::<f64>() - 0.5,
            rng.gen::<f64>() - 0.5,
        ) * 0.1; // Jitter amount to simulate area light

        let jittered_light_pos = light.position + light_jitter;
        let light_dir = (jittered_light_pos - ray.origin).unit_vector();
        let shadow_ray = Ray::new(ray.origin, light_dir);

        for sphere in spheres {
            if let Some(_) = sphere.hit(&shadow_ray) {
                return true;  // If any shadow ray hits an object, it's in shadow
            }
        }
    }
    false
}


fn random_in_unit_sphere() -> Vec3 {
    let mut rng = rand::thread_rng();
    loop {
        let p = Vec3::new(rng.gen::<f64>(), rng.gen::<f64>(), rng.gen::<f64>()) * 2.0 - Vec3::new(1.0, 1.0, 1.0);
        if p.dot(p) < 1.0 {
            return p;
        }
    }
}

fn color(ray: &Ray, spheres: &[Sphere], planes: &[Plane], light: &Light, depth: u32) -> Vec3 {
    if depth >= 50 {
        return Vec3::new(0.0, 0.0, 0.0); // Max recursion depth to avoid infinite loops
    }

    let mut closest_t = f64::MAX;
    let mut hit_color = None;

    // Check for sphere intersections
    for sphere in spheres {
        if let Some(t) = sphere.hit(ray) {
            if t < closest_t {
                closest_t = t;
                let hit_point = ray.point_at_parameter(t);
                let N = (hit_point - sphere.center).unit_vector();
                let light_dir = (light.position - hit_point).unit_vector();

                // Calculate new scattered ray for global illumination (random bounce)
                let target = hit_point + N + random_in_unit_sphere();
                let scattered_ray = Ray::new(hit_point, target - hit_point);
                let bounce_color = color(&scattered_ray, spheres, planes, light, depth + 1);

                // Combine direct lighting and bounce light (indirect lighting)
                hit_color = Some(bounce_color * 0.5 + Vec3::new(1.0, 0.5, 0.5) * N.dot(light_dir).max(0.0) * light.intensity);
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
                if in_shadow(&shadow_ray, spheres, &light) {
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
    // Sky color (if no object is hit)
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
    let aspect_ratio = nx as f64 / ny as f64;

    // Parameters for the camera
    let fov = 90.0;
    let aperture = 0.1;
    let focus_dist = 1.0;

    // Create a camera with the required parameters
    let camera = Camera::new(fov, aspect_ratio, aperture, focus_dist);

    // Create a light source
    let light = Light::new(Vec3::new(5.0, 5.0, 2.0), 1.0); // Point light

    // Create some spheres for rendering
    let spheres = vec![
        Sphere::new(Vec3::new(0.0, 0.0, -1.0), 0.5),
        Sphere::new(Vec3::new(1.0, 0.0, -1.0), 0.5),
        Sphere::new(Vec3::new(-1.0, 0.0, -1.0), 0.5),
    ];

    // Create a ground plane
    let planes = vec![
        Plane::new(Vec3::new(0.0, -0.5, 0.0), Vec3::new(0.0, 1.0, 0.0)), // Ground plane
    ];

    let mut img = RgbImage::new(nx as u32, ny as u32);
    let mut rng = rand::thread_rng();

    // Render the image
    for j in 0..ny {
        for i in 0..nx {
            let mut col = Vec3::new(0.0, 0.0, 0.0);
            for _ in 0..ns {
                let u = (i as f64 + rng.gen::<f64>()) / nx as f64;
                let v = (j as f64 + rng.gen::<f64>()) / ny as f64;
                let ray = camera.get_ray(u, v);
                col = col + color(&ray, &spheres, &planes, &light, 0);
            }
            col = col / ns as f64;

            // Convert color to 8-bit and store in the image buffer
            let ir = (255.99 * col.x) as u8;
            let ig = (255.99 * col.y) as u8;
            let ib = (255.99 * col.z) as u8;
            img.put_pixel(i as u32, (ny - j - 1) as u32, Rgb([ir, ig, ib]));
        }
    }

    // Save the image as PNG
    img.save("output.png").unwrap();
}
