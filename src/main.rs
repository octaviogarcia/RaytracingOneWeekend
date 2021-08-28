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

fn main() {
    //IMAGE
    let aspect_ratio: f64 = 16.0 / 9.0;
    let image_width:  u64  = 400;
    let image_height: u64 = ((image_width as f64) / aspect_ratio) as u64;
    
    //"Normal" world camera
    let camera = Camera::world_camera(90.,aspect_ratio);
    //Weird camera
    //let camera = Camera::new(Point3::new(-2.,2.,1.),Point3::new(0.,0.,-1.),Vec3::new(0.,1.,0.),90.,aspect_ratio,0.,1.);
    //let camera = Camera::new(Point3::new(-2.,2.,1.),Point3::new(0.,0.,-1.),Vec3::new(0.,1.,0.),20.,aspect_ratio,0.,1.);
    //12.2 
    let lookfrom = Point3::new(3.,3.,2.);
    let lookat   = Point3::new(0.,0.,-1.);
    let vup      = Vec3::new(0.,1.,0.);
    let aperture = 2.;
    let dist_to_focus = (lookfrom-lookat).length();
    let camera = Camera::new(lookfrom,lookat,vup,20.,aspect_ratio,aperture,dist_to_focus);

    let samples_per_pixel = 100;
    let max_depth = 50;

    //OBJECTS
    let mut world = HittableList::new();
    {//10.5
        //Creating objects with only 1 instance is bad practice, should be a function really...
        let mat_ground = Rc::new(Lambertian::new(Color::new(0.8,0.8,0.0)));
        let mat_center = Rc::new(Lambertian::new(Color::new(0.1,0.2,0.5)));
        let mat_left   = Rc::new(Dieletric::new(1.5));
        let mat_right  = Rc::new(Metal::new_fuzz(Color::new(0.8,0.6,0.2),0.0));
        world.add(Box::new(Sphere{center: Point3::new( 0.0, -100.5, -1.0), radius: 100.0, material: mat_ground}));
        world.add(Box::new(Sphere{center: Point3::new( 0.0,    0.0, -1.0), radius:   0.5, material: mat_center}));
        world.add(Box::new(Sphere{center: Point3::new(-1.0,    0.0, -1.0), radius:   0.5, material: mat_left.clone()}));
        //10.5
        //world.add(Box::new(Sphere{center: Point3::new(-1.0,    0.0, -1.0), radius:  -0.4, material: mat_left.clone()}));
        //11.2
        world.add(Box::new(Sphere{center: Point3::new(-1.0,    0.0, -1.0), radius:  -0.45, material: mat_left.clone()}));
        world.add(Box::new(Sphere{center: Point3::new( 1.0,    0.0, -1.0), radius:   0.5, material: mat_right}));
    }
    /*{//11.1
        let mat_left  = Rc::new(Lambertian::new(Color::new(0.,0.,1.)));
        let mat_right = Rc::new(Lambertian::new(Color::new(1.,0.,0.)));
        let r = (PI/4.0).cos();
        world.add(Box::new(Sphere{center: Point3::new(-r, 0.,-1.0), radius: r, material: mat_left}));
        world.add(Box::new(Sphere{center: Point3::new( r, 0.,-1.0), radius: r, material: mat_right}));
    }*/

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
