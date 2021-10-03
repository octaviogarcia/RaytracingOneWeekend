extern crate num_cpus;
extern crate sdl2;

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
use materials::*;

use sdl2::event::Event;
use sdl2::keyboard::Keycode;
use std::time::Duration;
use rand::prelude::*;

#[derive(Copy,Clone)]
struct ColorsBox {
    pub colors: *mut Vec<Color>
}
unsafe impl Send for ColorsBox{}

fn ray_color(r: &Ray,world: &HittableList, depth: u64,tmin: f64,tmax: f64) -> Color{
    if depth == 0 {
        return Color::ZERO;
    }
    match world.hit(r,tmin,tmax) {
        Some(hr) => {
            let rslt = hr.material.scatter(r,&hr);
            return rslt.attenuation*ray_color(&rslt.ray, world,depth-1,tmin,tmax);
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
    let mat_ground = Material::new_lambertian(Color::new(0.5,0.5,0.5));
    //world.add_marched_sphere(&MarchedSphere{center: Point3::new(0., -1000.,0.), radius: 1000.0, material: mat_ground});
    world.add_sphere(&Sphere{center: Point3::new(0., -1000.,0.), radius: 1000.0, material: mat_ground});
    for a in -11..11{
        let af = a as f64;
        for b in -11..11{
            let bf = b as f64;
            let center = Point3::new(af+0.9*f64::rand(),0.2,bf+0.9*f64::rand());
            let add_to_world = (center - Point3::new(4.,0.2,0.)).length() > 0.9;
            if add_to_world{
                let sphere_mat: Material;
                let mat_prob = f64::rand();
                if mat_prob < 0.8{
                    let albedo = Color::rand() * Color::rand();
                    sphere_mat = Material::new_lambertian(albedo);
                }
                else if mat_prob < 0.95{
                    let albedo = Color::rand_range(0.5,1.);
                    let fuzz   = f64::rand_range(0.,0.5);
                    sphere_mat = Material::new_metal_fuzz(albedo,fuzz);
                }
                else{
                    sphere_mat = Material::new_dielectric(1.5);
                }
                world.add_sphere(&Sphere{center: center, radius: 0.2, material: sphere_mat});
            }
        }
    }
    {
        let mat = Material::new_dielectric(1.5);
        //world.add_sphere(&Sphere{center: Point3::new(0.,1.,0.), radius: 1., material: mat});
        //world.add_marched_sphere(&MarchedSphere{center: Point3::new(0.,1.,0.), radius: 1., material: mat});
        world.add_marched_box(&MarchedBox{center: Point3::new(0.,1.,0.), sizes: Vec3::new(0.5,0.5,0.5), material: mat});
    }
    {
        let mat = Material::new_lambertian(Color::new(0.4,0.2,0.1));
        //world.add_sphere(&Sphere{center: Point3::new(-4.,1.,0.), radius: 1., material: mat});
        //world.add_marched_sphere(&MarchedSphere{center: Point3::new(-4.,1.,0.), radius: 1., material: mat});
        world.add_marched_box(&MarchedBox{center: Point3::new(-4.,1.,0.), sizes: Vec3::new(0.3,0.3,0.3), material: mat});
    }
    {
        let mat = Material::new_metal(Color::new(0.7,0.6,0.5));
        //world.add_sphere(&Sphere{center: Point3::new(4.,1.,0.), radius: 1., material: mat});
        //world.add_marched_sphere(&MarchedSphere{center: Point3::new(4.,1.,0.), radius: 1., material: mat});
        world.add_marched_box(&MarchedBox{center: Point3::new(4.,1.,0.), sizes: Vec3::new(0.5,0.5,0.5), material: mat});
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

use std::thread;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

fn draw(camera: &Camera,world: &HittableList,max_depth: u64,tmin: f64,tmax: f64,
    samples_per_pixel: u64,image_width: u64,image_height: u64,
    colors_box: ColorsBox,
    tid: u64,assigned_thread: &Vec<u64>,pixels: &AtomicU64)
{
    let image_width_f  = image_width as f64;
    let image_height_f = image_height as f64;

    for line in 0..image_height {
        let j_f = line as f64;
        for col in 0..image_width{
            let pos = (line*image_width+col) as usize;
            if assigned_thread[pos] != tid {continue;}

            let i_f = col as f64;
            let mut pixel_color = Color::ZERO;
            for _s in 0..samples_per_pixel{
                let u = (i_f+f64::rand())/(image_width_f-1.);
                let v = (j_f+f64::rand())/(image_height_f-1.);
                let ray = camera.get_ray(u,1.0-v);
                pixel_color += ray_color(&ray,&world,max_depth,tmin,tmax);
            }
            unsafe { (*(colors_box.colors))[pos] = normalize_color(pixel_color,samples_per_pixel); }

            pixels.fetch_add(1,Ordering::Relaxed);//Inform pixel is done to Log Thread
        }
    }
}

fn main() {
    //IMAGE
    let aspect_ratio: f64 = 3.0 / 2.0;
    let image_width:    u64 = 1000;
    let image_width_f:  f64 = image_width as f64;
    let image_height_f: f64 = image_width_f/ aspect_ratio;
    let image_height:   u64 = image_height_f as u64;
    let image_size:     u64 = image_width*image_height;
    let image_size_f:   f64 = image_width_f*image_height_f;

    let camera: Camera;
    {
        let lookfrom = Point3::new(13.,2.,3.);
        let lookat   = Point3::new(0.,0.,0.);
        let vup      =   Vec3::new(0.,1.,0.);
        let aperture = 0.1;
        let dist_to_focus = 10.;
        camera = Camera::new(lookfrom,lookat,vup,20.,aspect_ratio,aperture,dist_to_focus);
    }

    let samples_per_pixel = 20;
    let max_depth = 20;
    let world = random_scene();

    let pixels_atomic = AtomicU64::new(0);
    let arc_pixels_atomic = Arc::new(pixels_atomic);
    
    {//Log thread
        let pxls_atom = arc_pixels_atomic.clone();
        thread::spawn(move || {
            loop {
                let progress = pxls_atom.load(Ordering::Relaxed);
                print_progress((progress as f64)/image_size_f);
                if image_size == progress{ 
                    print_progress(1.0);
                    return;
                }
                ::std::thread::sleep(Duration::new(0, 1_000_000_000u32 / 2));
            }
        });
    }

    let num_threads = (num_cpus::get()-1) as u64;//Leave 1 thread for logging and drawing

    let mut assigned_thread: Vec<u64> = Vec::with_capacity(image_size as usize);
    for cidx in 0..image_size{
        assigned_thread.push(cidx%num_threads);
    }
    assigned_thread.shuffle(&mut rand::thread_rng());
    assigned_thread.shrink_to_fit();
    let arc_assigned_thread = Arc::new(assigned_thread);

    let colors_box = ColorsBox{colors: &mut vec!(-Color::ZERO;image_size as usize)};

    let mut handlers: Vec<thread::JoinHandle<()>> = Vec::with_capacity(num_threads as usize);
    let arc_camera = Arc::new(camera);
    let arc_world = Arc::new(world);
    eprintln!("Running {} threads",num_threads);

    for i in 0..num_threads {
        let cam = arc_camera.clone();
        let wrld = arc_world.clone();
        let pxls_atom = arc_pixels_atomic.clone();
        let assgn_th = arc_assigned_thread.clone();
        let tmin = 0.001;
        let tmax = 100.0;//@TODO: You could find these from bounding boxes from the scene
        let draw_thread = move || {
            return draw(&cam,&wrld,max_depth,tmin,tmax,
                samples_per_pixel,image_width,image_height,
                colors_box,
                i,&assgn_th,&pxls_atom);
        };
        handlers.push(thread::spawn(draw_thread));
    }

    let colors: &mut Vec<Color> = unsafe {&mut (*colors_box.colors) };
    draw_to_sdl(&colors,image_width,image_height);
    /*
    {
        for h in handlers{
            h.join().unwrap();
        }
        arc_pixels_atomic.clone().store(image_size,Ordering::Relaxed);
        write_ppm(&colors,image_width,image_height);
    }*/
}

fn draw_to_sdl(colors: &Vec<Color>,image_width: u64,image_height: u64){
    let sdl_context = sdl2::init().unwrap();
    let video_subsystem = sdl_context.video().unwrap();

    let window = video_subsystem.window("raytracer", image_width as u32, image_height as u32)
    .position_centered().build().unwrap();

    let mut canvas = window.into_canvas().build().unwrap();
    canvas.set_draw_color(sdl2::pixels::Color::RGB(0,0,0));
    canvas.clear();
    canvas.present();

    let mut event_pump = sdl_context.event_pump().unwrap();
    let mut pixels_run = 0;//Cache the longest pixel run... to avoid retexturing too much. Requires initialization to -0.0
    let texture_creator = canvas.texture_creator();
    let mut texture = texture_creator
    .create_texture_target(texture_creator.default_pixel_format(), image_width as u32, image_height as u32)
    .unwrap();

    'running: loop {
        canvas.with_texture_canvas(&mut texture, |texture_canvas| {
            let mut run = true;
            for pos in pixels_run..(image_width*image_height){
                let c = colors[pos as usize];

                let is_zero = c.x() == 0.;//@HACK: Expects intialization to -0.0, so we know if it actually was rendered or is just 0.
                let is_neg = (1. as f64).copysign(c.x()) == -1.;
                run = run && !(is_zero && is_neg);
                pixels_run += run as u64;

                texture_canvas.set_draw_color(sdl2::pixels::Color::RGB((c.x()*256.0) as u8,(c.y()*256.0) as u8,(c.z()*256.0) as u8));
                let y = pos / image_width;
                let x = pos - y*image_width;
                texture_canvas.draw_point(sdl2::rect::Point::new(x as i32, y as i32)).unwrap();
            }
        }).unwrap();

        canvas.clear();
        for event in event_pump.poll_iter() {
            match event {
                Event::Quit {..} |
                Event::KeyDown { keycode: Some(Keycode::Escape), .. } => {
                    break 'running
                },
                _ => {}
            }
        }
        canvas.copy(&texture,None,None).unwrap();
        canvas.present();
        ::std::thread::sleep(Duration::new(0, 1_000_000_000u32 / 15));
    }
}