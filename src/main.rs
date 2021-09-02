extern crate num_cpus;

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
use materials::Dieletric;

fn ray_color(r: &Ray,world: &HittableList, depth: u64) -> Color{
    if depth == 0 {
        return Color::ZERO;
    }
    match world.hit(r,0.001,INF) {
        Some(hr) => {
            let rslt = hr.material.scatter(r,&hr);
            return rslt.attenuation*ray_color(&rslt.ray, world,depth-1);
        },
        None => {
            let unit_dir: Vec3 = r.dir.unit();
            let t: f64 = 0.5*(unit_dir.y() + 1.0);
            return lerp(t,Color::new(1.0,1.0,1.0),Color::new(0.5,0.7,1.0));
        },
    }
}


fn random_scene() -> HittableList{
    let mut world = HittableList::new();
    let mat_ground = Arc::new(Lambertian::new(Color::new(0.5,0.5,0.5)));
    world.add_sphere(Box::new(Sphere{center: Point3::new(0., -1000.,0.), radius: 1000.0, material: mat_ground}));
    for a in -11..11{
        let af = a as f64;
        for b in -11..11{
            let bf = b as f64;
            let center = Point3::new(af+0.9*f64::rand(),0.2,bf+0.9*f64::rand());
            let add_to_world = (center - Point3::new(4.,0.2,0.)).length() > 0.9;
            if add_to_world{
                let sphere_mat: Arc<dyn Material + Send + Sync>;
                let mat_prob = f64::rand();
                if mat_prob < 0.8{
                    let albedo = Color::rand() * Color::rand();
                    sphere_mat = Arc::new(Lambertian::new(albedo));
                }
                else if mat_prob < 0.95{
                    let albedo = Color::rand_range(0.5,1.);
                    let fuzz   = f64::rand_range(0.,0.5);
                    sphere_mat = Arc::new(Metal::new_fuzz(albedo,fuzz));
                }
                else{
                    sphere_mat = Arc::new(Dieletric::new(1.5));
                }
                world.add_sphere(Box::new(Sphere{center: center, radius: 0.2, material: sphere_mat}));
            }
        }
    }
    {
        let mat = Arc::new(Dieletric::new(1.5));
        world.add_sphere(Box::new(Sphere{center: Point3::new(0.,1.,0.), radius: 1., material: mat}));
    }
    {
        let mat = Arc::new(Lambertian::new(Color::new(0.4,0.2,0.1)));
        world.add_sphere(Box::new(Sphere{center: Point3::new(-4.,1.,0.), radius: 1., material: mat}));
    }
    {
        let mat = Arc::new(Metal::new(Color::new(0.7,0.6,0.5)));
        world.add_sphere(Box::new(Sphere{center: Point3::new(4.,1.,0.), radius: 1., material: mat}));
    }
    return world;
}

fn round_n(f: f64,n: u64) -> f64 {
    let mut mul = 1.0;
    for _n in 0..n{
        mul *= 10.0;
    }
    return (f*mul).round()/mul;
}



fn print_progress(progress: f64) -> (){
    let progress100 = round_n(100.0*progress,2);
    let frac = progress100 % 1.;
    let int =  progress100 - frac;
    eprint!("{:>3}.{:0>2}%\r",(int as u64),((frac*100.) as u64));
}

use std::{thread,time};
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

fn draw(camera: &Camera,world: &HittableList,max_depth: u64,samples_per_pixel: u64,image_width: u64,image_height: u64,samples: &AtomicU64) -> Vec<Color>{
    let image_width_f  = image_width as f64;
    let image_height_f = image_height as f64;
    let mut colors: Vec<Color> = vec!(Color::ZERO;(image_height*image_width) as usize);
    for line in 0..image_height {
        let j_f = line as f64;
        for col in 0..image_width{
            let i_f = col as f64;
            let mut pixel_color = Color::ZERO;
            for _s in 0..samples_per_pixel{
                let u = (i_f+f64::rand())/(image_width_f-1.);
                let v = (j_f+f64::rand())/(image_height_f-1.);
                let ray = camera.get_ray(u,1.0-v);
                pixel_color += ray_color(&ray,&world,max_depth);
            }
            samples.fetch_add(samples_per_pixel,Ordering::Relaxed);
            let pos = line*image_width+col;
            colors[pos as usize] = pixel_color;
        }
    }
    return colors;
}


fn main() {
    //IMAGE
    let aspect_ratio: f64 = 3.0 / 2.0;
    let image_width:    u64 = 400;
    let image_width_f:  f64 = image_width as f64;
    let image_height_f: f64 = image_width_f/ aspect_ratio;
    let image_height:   u64 = image_height_f as u64;

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
    let max_depth = 20;
    let world = random_scene();

    let samples_atomic = AtomicU64::new(0);
    let arc_samples_atomic = Arc::new(samples_atomic);
    let total_samples = image_height*image_width*samples_per_pixel;
    
    {//Log thread
        let smpls_atomic = arc_samples_atomic.clone();
        thread::spawn(move || {
            loop {
                let progress = smpls_atomic.load(Ordering::Relaxed);
                print_progress((progress as f64)/(total_samples as f64));
                if total_samples == progress{ 
                    return;
                }
                thread::sleep(time::Duration::from_millis(100));
            }
        });
    }

    let num_threads = num_cpus::get() as u64;
    let mut handlers: Vec<thread::JoinHandle<Vec<Color>>> = Vec::with_capacity(num_threads as usize);
    let arc_camera = Arc::new(camera);
    let arc_world = Arc::new(world);
    let samples_per_pixel_per_thread = samples_per_pixel / num_threads;
    let missing_samples_per_pixel    = samples_per_pixel % num_threads;
    eprintln!("Running {} threads",num_threads);
    for i in 0..num_threads {
        let cam = arc_camera.clone();
        let wrld = arc_world.clone();
        let smpls_atomic = arc_samples_atomic.clone();

        let samples: u64;
        //Load the starting threads with an extra pixel so we get exactly what se set at the top
        if i < missing_samples_per_pixel {
            samples = samples_per_pixel_per_thread + 1;
        }
        else{
            samples = samples_per_pixel_per_thread;
        }

        let draw_thread = move || {
            return draw(&cam,&wrld,max_depth,samples,image_width,image_height,&smpls_atomic);
        };
        handlers.push(thread::spawn(draw_thread));
    }

    let mut colors: Vec<Color> = vec!(Color::ZERO;(image_height*image_width) as usize);
    for h in handlers{
        //@Speed: This blocks threads sequentially from the order initialization, this order of join() is not optimal
        //Rust has no try_join, so we need to work around it with channels... maybe once this is simpler I will complicate it
        let local_clrs = h.join().unwrap();
        for cidx in 0..local_clrs.len(){
            colors[cidx] += local_clrs[cidx];
        }
    }

    for cidx in 0..colors.len(){
        colors[cidx] = normalize_color(colors[cidx],samples_per_pixel);
    }

    write_ppm(&colors,image_width,image_height);
    arc_samples_atomic.clone().store(total_samples,Ordering::Relaxed);
}
