mod vec3;
use vec3::*;

mod utils;
use utils::*;

mod ray;
use ray::*;

mod hits;
use hits::*;

mod camera;
use camera::*;


fn ray_color(r: &Ray,world: &HittableList) -> Color{
    match world.hit(r,0.,INF) {
        None => {
            let unit_dir: Vec3 = r.dir.unit();
            let t: f64 = 0.5*(unit_dir.y() + 1.0);
            return lerp(t,Color::new(1.0,1.0,1.0),Color::new(0.5,0.7,1.0));
        },
        Some(hr) => {
            let offset = Color::new(0.5,0.5,0.5);
            return 0.5*hr.normal + offset;
        }
    }
}

fn main() {
    //IMAGE
    let aspect_ratio: f64 = 16.0 / 9.0;
    let image_width:  u64  = 400;
    let image_height: u64 = ((image_width as f64) / aspect_ratio) as u64;
    
    //VIEWPORT
    let viewport_height: f64 = 2.0;
    let viewport_width:  f64 = viewport_height * aspect_ratio;
    let focal_length:    f64 = 1.0;

    let camera = Camera::new(viewport_height,viewport_width,focal_length);
    let samples_per_pixel = 10;

    //OBJECTS
    let mut world = HittableList::new();
    world.add(Box::new(Sphere{center: Point3::new(0.,0.,-1.),radius: 0.5}));
    world.add(Box::new(Sphere{center: Point3::new(0.,-100.5,-1.),radius: 100.}));

    eprintln!("Inicio");
    println!("P3\n{} {}\n255",image_width,image_height);
    for j in (0..image_height).rev() {
        eprint!("Faltan {} lineas        \r",j);
        for i in 0..image_width{
            let mut pixel_color = Color::ZERO;
            for _s in 0..samples_per_pixel{
                let u = ((i as f64)+rand_f64())/(image_width as f64-1.);
                let v = ((j as f64)+rand_f64())/(image_height as f64-1.);
                let ray = camera.get_ray(u,v);
                pixel_color += ray_color(&ray,&world);
            }
            write_color(&pixel_color,samples_per_pixel);
        }
    }
    eprintln!("                         ");
    eprintln!("Finalizo");
}
