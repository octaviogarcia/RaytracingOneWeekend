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

mod materials;
use materials::Material;
use materials::Lambertian;
use materials::Metal;
use materials::MaterialScatterResult;
use materials::Dieletric;

fn ray_color(r: &Ray,world: &HittableList, depth: u64) -> Color{
    if depth == 0 {
        return Color::ZERO;
    }
    match world.hit(r,0.001,INF) {
        Some(hr) => {
            let rslt = hr.material.scatter(r,&hr);
            match rslt {
                Some(r) => {
                    return r.attenuation*ray_color(&r.ray, world,depth-1);
                }
                None => { return Color::ZERO; } //Nevear reached by Lambertian, Metal
            }
        },
        None => {
            let unit_dir: Vec3 = r.dir.unit();
            let t: f64 = 0.5*(unit_dir.y() + 1.0);
            return lerp(t,Color::new(1.0,1.0,1.0),Color::new(0.5,0.7,1.0));
        },
    }
}

use std::rc::Rc;

fn random_scene() -> HittableList{
    let mut world = HittableList::new();
    let mat_ground = Rc::new(Lambertian::new(Color::new(0.5,0.5,0.5)));
    world.add(Box::new(Sphere{center: Point3::new(0., -1000.,0.), radius: 1000.0, material: mat_ground}));
    for a in -11..11{
        let af = a as f64;
        for b in -11..11{
            let bf = b as f64;
            let center = Point3::new(af+0.9*f64::rand(),0.2,bf+0.9*f64::rand());
            let add_to_world = (center - Point3::new(4.,0.2,0.)).length() > 0.9;
            if add_to_world{
                let sphere_mat: Rc<dyn Material>;
                let mat_prob = f64::rand();
                if mat_prob < 0.8{
                    let albedo = Color::rand() * Color::rand();
                    sphere_mat = Rc::new(Lambertian::new(albedo));
                }
                else if mat_prob < 0.95{
                    let albedo = Color::rand_range(0.5,1.);
                    let fuzz   = f64::rand_range(0.,0.5);
                    sphere_mat = Rc::new(Metal::new_fuzz(albedo,fuzz));
                }
                else{
                    sphere_mat = Rc::new(Dieletric::new(1.5));
                }
                world.add(Box::new(Sphere{center: center, radius: 0.2, material: sphere_mat}));
            }
        }
    }
    {
        let mat = Rc::new(Dieletric::new(1.5));
        world.add(Box::new(Sphere{center: Point3::new(0.,1.,0.), radius: 1., material: mat}));
    }
    {
        let mat = Rc::new(Lambertian::new(Color::new(0.4,0.2,0.1)));
        world.add(Box::new(Sphere{center: Point3::new(-4.,1.,0.), radius: 1., material: mat}));
    }
    {
        let mat = Rc::new(Metal::new(Color::new(0.7,0.6,0.5)));
        world.add(Box::new(Sphere{center: Point3::new(4.,1.,0.), radius: 1., material: mat}));
    }
    return world;
}

fn main() {
    //IMAGE
    let aspect_ratio: f64 = 3.0 / 2.0;
    let image_width:  u64  = 1200;
    let image_height: u64 = ((image_width as f64) / aspect_ratio) as u64;

    let camera: Camera;
    {
        let lookfrom = Point3::new(13.,2.,3.);
        let lookat   = Point3::new(0.,0.,0.);
        let vup      =   Vec3::new(0.,1.,0.);
        let aperture = 0.1;
        let dist_to_focus = 10.;
        camera = Camera::new(lookfrom,lookat,vup,20.,aspect_ratio,aperture,dist_to_focus);
    }

    let samples_per_pixel = 500;
    let max_depth = 50;
    
    let world = random_scene();

    eprintln!("Inicio");
    println!("P3\n{} {}\n255",image_width,image_height);
    for j in (0..image_height).rev() {
        eprint!("Faltan {} lineas        \r",j);
        for i in 0..image_width{
            let mut pixel_color = Color::ZERO;
            for _s in 0..samples_per_pixel{
                let u = ((i as f64)+f64::rand())/(image_width as f64-1.);
                let v = ((j as f64)+f64::rand())/(image_height as f64-1.);
                let ray = camera.get_ray(u,v);
                pixel_color += ray_color(&ray,&world,max_depth);
            }
            write_color(&pixel_color,samples_per_pixel);
        }
    }
    eprintln!("                         ");
    eprintln!("Finalizo");
}
