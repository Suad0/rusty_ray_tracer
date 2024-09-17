mod bsp;
mod ray_tracer;


fn main() {
    /*
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


     */

    let mut tree = crate::bsp::BSPTree::new();
    tree.add_point(crate::bsp::Point { x: 5.0, y: 3.0 });
    tree.add_point(crate::bsp::Point { x: 1.0, y: 7.0 });
    tree.add_point(crate::bsp::Point { x: 10.0, y: 2.0 });
    tree.add_point(crate::bsp::Point { x: 3.0, y: 6.0 });
    tree.add_point(crate::bsp::Point { x: 8.0, y: 1.0 });

    println!("BSP Tree points:");
    tree.print();
}
