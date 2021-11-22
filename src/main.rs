extern crate num_cpus;
extern crate sdl2;

mod math;
use math::vec3::*;

mod utils;
use utils::*;

mod ray;

mod hits;
use hits::*;

mod traced;
use traced::*;

mod marched;
use marched::*;

mod camera;
use camera::*;

mod materials;
use materials::*;

mod render_thread;

use sdl2::event::Event;
use sdl2::keyboard::Keycode;
use std::time::Duration;
use crate::math::mat4x4::Mat4x4;

fn random_scene() -> HittableList{
    let mut world = HittableList::new();
    let mat_ground = Material::new_lambertian(Color::new(0.5,0.5,0.5));
    //world+=&MarchedSphere{center: Point3::new(0., -1000.,0.), radius: 1000.0, material: mat_ground};
    world+=&Sphere::new_with_radius(&Point3::new(0., -1000.,0.),1000.0,&mat_ground);
    for a in -11..11{
        let af = a as f32;
        for b in -11..11{
            let bf = b as f32;
            let center = Point3::new(af+0.9*f32::rand(),0.2,bf+0.9*f32::rand());
            let add_to_world = (center - Point3::new(4.,0.2,0.)).length() > 0.9;
            if add_to_world{
                let sphere_mat: Material;
                let mat_prob = f32::rand();
                if mat_prob < 0.8{
                    let albedo = Color::rand() * Color::rand();
                    sphere_mat = Material::new_lambertian(albedo);
                }
                else if mat_prob < 0.95{
                    let albedo = Color::rand_range(0.5,1.);
                    let fuzz   = f32::rand_range(0.,0.5);
                    sphere_mat = Material::new_metal_fuzz(albedo,fuzz);
                }
                else{
                    sphere_mat = Material::new_dielectric(1.5);
                }

                //READ BOTTOM UP in order of operations
                let m_local_to_world = m4x4!(TR center)
                ^m4x4!(RX f32::rand()*2.*PI)^m4x4!(RY f32::rand()*2.*PI)^m4x4!(RZ f32::rand()*2.*PI)//Rotate randomly
                ^m4x4!(SC f32::rand()+1.,f32::rand()+1.,f32::rand()+1.)//Warp into an egg
                ^m4x4!(SC 0.2,0.2,0.2);//Set Radius
                world+=&Sphere::new(&m_local_to_world,&sphere_mat);
            }
        }
    }
    {
        let mat = Material::new_dielectric(1.5);
        //world+=&Sphere{center: Point3::new(0.,1.,0.), radius: 1., material: mat};
        //world+=&MarchedSphere{center: Point3::new(0.,1.,0.), radius: 1., material: mat};
        //world+=&MarchedBox{center: Point3::new(0.,1.,0.), sizes: Vec3::new(0.5,0.5,0.5), material: mat};
        //world+=Arc::new(MarchedTorus{center: Point3::new(0.,1.,0.), sizes: Vec3::new(0.5,0.1,0.1), material: mat}) as Arc<dyn Marched + Send + Sync>;
        world+=&MarchedTorus{center: Point3::new(0.,1.,0.), sizes: Vec3::new(0.5,0.1,0.1), material: mat};
        //world+=&InfinitePlane{center: Point3::new(0.,1.,0.),normal: Vec3::new(0.,0.,1.), material: mat};
    }
    {
        //let mat = Material::new_lambertian(Color::new(0.4,0.2,0.1));
        //world+=&Sphere{center: Point3::new(-4.,1.,0.), radius: 1., material: mat};
        //world+=&MarchedSphere{center: Point3::new(-4.,1.,0.), radius: 1., material: mat};
        //world+=&MarchedBox{center: Point3::new(-4.,1.,0.), sizes: Vec3::new(0.3,0.3,0.3), material: mat};
        //world+=&InfinitePlane::new(&Point3::new(-4.,1.,0.),&Vec3::new(0.,0.,1.),&mat};
        let p1 = Point3::new(7.,1.,0.);
        let p2 = Point3::new(6.,1.1,0.5);
        let p3 = Point3::new(6.,1.5,0.);
        world+=&Parallelogram::new3points(
            &p1,&p2,&p3,&Material::new_metal(Color::new(1.,0.5,1.))
        );
        world+=&Triangle::new3points(
            &(p1+Point3::new(0.,0.5,0.)),&p2,&p3,&Material::new_lambertian(Color::new(1.,1.,0.))
        );
    }
    {
        let mat = Material::new_metal(Color::new(0.7,0.6,0.5));
        //world+=&Sphere{center: Point3::new(4.,1.,0.), radius: 1., material: mat};
        //world+=&MarchedSphere{center: Point3::new(4.,1.,0.), radius: 1., material: mat};
        //world+=&MarchedBox{center: Point3::new(4.,1.,0.), sizes: Vec3::new(0.5,0.5,0.5), material: mat};
        let m_local_to_world = m4x4!(TR 4.,1.,0.)//Move it
        ^m4x4!(RX f32::rand()*2.*PI)^m4x4!(RY f32::rand()*2.*PI)^m4x4!(RZ f32::rand()*2.*PI);//Rotate randomly
        world+=&Cube::new(&m_local_to_world,&mat);
        //world+=&Cube::new_with_length(&Point3::new(4.,1.,0.),1.,&mat);
        //world+=&InfinitePlane{center: Point3::new(4.,1.,0.),normal: Vec3::new(0.,0.,1.), material:  Material::new_metal(Color::new(1.,1.,1.))};
    }
    return world;
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

fn main() {
    //IMAGE
    let aspect_ratio: f32 = 3.0 / 2.0;
    let image_width:    u32 = 1000;
    let image_width_f:  f32 = image_width as f32;
    let image_height_f: f32 = image_width_f/ aspect_ratio;
    let image_height:   u32 = image_height_f as u32;
    let image_size:     u32 = image_width*image_height;

    let camera: Camera;
    {
        let lookfrom = Point3::new(13.,2.,3.);
        //let lookfrom = Point3::new(13.,1.,0.);
        let lookat   = Point3::new(0.,0.,0.);
        //let lookat   = Point3::new(0.,1.,0.);
        let vup      =   Vec3::new(0.,1.,0.);
        let aperture = 0.1;
        let dist_to_focus = 10.;
        camera = Camera::new(lookfrom,lookat,vup,20.,aspect_ratio,aperture,dist_to_focus);
    }

    let samples_per_pixel: u32 = 200;
    let max_depth: u32 = 50;
    let world = random_scene();

    let samples_atomic = AtomicU64::new(0);
    let arc_samples_atomic = Arc::new(samples_atomic);
    
    {//Log thread
        let smpls_atom = arc_samples_atomic.clone();
        thread::spawn(move || {
            let total_samples = (image_size as u64)*(samples_per_pixel as u64);
            let total_samples_f = total_samples as f64;
            let start = std::time::Instant::now();
            loop {
                let progress = smpls_atom.load(Ordering::Relaxed);
                print_progress((progress as f64)/total_samples_f);
                if total_samples == progress{ 
                    print_progress(1.0);
                    break;
                }
                ::std::thread::sleep(Duration::new(0, 1_000_000_000u32 / 2));
            }
            eprintln!("{} seconds",start.elapsed().as_secs());
        });
    }

    let num_threads = num_cpus::get() as u32 - 1;

    let mut assigned_thread: Vec<u32> = Vec::with_capacity(image_size as usize);
    {
        const CACHE_SIZE: u32 = 32*1024;
        const CHUNK_SIZE: u32 = CACHE_SIZE/std::mem::size_of::<Color>() as u32;
        for chunk in 0..(image_size / CHUNK_SIZE){
            let id = chunk % num_threads;
            for _i in 0..CHUNK_SIZE{
                assigned_thread.push(id);
            }
        }
        {
            let id = (image_size / CHUNK_SIZE) % num_threads;//Last iteration is (image_size / CHUNK_SIZE - 1) % num_threads
            let leftover = image_size - image_size/CHUNK_SIZE*CHUNK_SIZE;
            for _i in 0..leftover{
                assigned_thread.push(id);
            }
        }
    }
    assigned_thread.shrink_to_fit();
    let arc_assigned_thread = Arc::new(assigned_thread);

    let colors_box = render_thread::ColorsBox{colors: &mut vec!(Color::ZERO;image_size as usize),
                               samples: &mut vec!(0;image_size as usize),
                               true_samples: &mut vec!(0;image_size as usize)};

    let mut handlers: Vec<thread::JoinHandle<()>> = Vec::with_capacity(num_threads as usize);
    let arc_camera = Arc::new(camera);
    let arc_world = Arc::new(world);
    eprintln!("Running {} threads",num_threads);

    for i in 0..num_threads {
        let cam = arc_camera.clone();
        let wrld = arc_world.clone();
        let smpls_atom = arc_samples_atomic.clone();
        let assgn_th = arc_assigned_thread.clone();
        let tmin = 0.001;
        let tmax = 100.0;//@TODO: You could find these from bounding boxes from the scene
        
        let draw_thread = move || {
            return render_thread::render(&cam,&wrld.freeze(&cam),max_depth,tmin,tmax,
                samples_per_pixel,image_width,image_height,
                colors_box,
                i,&assgn_th,&smpls_atom);
        };
        handlers.push(thread::spawn(draw_thread));
    }

    let colors: &mut Vec<Color> = unsafe {&mut (*colors_box.colors) };
    let samples: &mut Vec<u32> = unsafe {&mut (*colors_box.samples) };
    let true_samples: &mut Vec<u32> = unsafe {&mut (*colors_box.true_samples) };
    draw_to_sdl(&colors,&samples,&true_samples,samples_per_pixel,image_width,image_height);
    /*
    {
        for h in handlers{
            h.join().unwrap();
        }
        arc_pixels_atomic.clone().store(image_size,Ordering::Relaxed);
        write_ppm(&colors,samples_per_pixel,image_width,image_height);
    }*/
}

fn draw_to_sdl(colors: &Vec<Color>,samples: &Vec<u32>,true_samples: &Vec<u32>,_samples_per_pixel: u32,image_width: u32,image_height: u32){
    let sdl_context = sdl2::init().unwrap();
    let video_subsystem = sdl_context.video().unwrap();

    let window = video_subsystem.window("ottomarcher", image_width as u32, image_height as u32)
    .position_centered().build().unwrap();

    let mut canvas = window.into_canvas().build().unwrap();
    canvas.set_draw_color(sdl2::pixels::Color::RGB(0,0,0));
    canvas.clear();
    canvas.present();

    let mut event_pump = sdl_context.event_pump().unwrap();
    let texture_creator = canvas.texture_creator();
    let mut show_samples = false;
    let mut pixels = vec!(0 as u8;(image_width*image_height*3) as usize);
    'running: loop {
        if show_samples{
            let mut max_samples = 1.0;
            for pos in 0..image_width*image_height{
                if (true_samples[pos as usize] as f32) > max_samples {
                    max_samples = true_samples[pos as usize] as f32;
                }
            }
            for pos in 0..image_width*image_height{
                let aux = (pos*3) as usize;
                let smpls = (true_samples[pos as usize] as f32)/max_samples;
                let c = normalize_color(Color::new(smpls,smpls,smpls),1);
                pixels[aux+0] = (c.x()*256.0) as u8;
                pixels[aux+1] = (c.y()*256.0) as u8;
                pixels[aux+2] = (c.z()*256.0) as u8;
            }
        }
        else{
            for pos in 0..image_width*image_height{
                let aux = (pos*3) as usize;
                let smpls = samples[pos as usize];
                let c = normalize_color(colors[pos as usize],smpls);
                pixels[aux+0] = (c.x()*256.0) as u8;
                pixels[aux+1] = (c.y()*256.0) as u8;
                pixels[aux+2] = (c.z()*256.0) as u8;
            }
        }
        //pitch = row in bytes. 1 byte per color -> 3*width
        let texture = sdl2::surface::Surface::from_data(pixels.as_mut_slice(), image_width, image_height, image_width*3, sdl2::pixels::PixelFormatEnum::RGB24)
        .unwrap().as_texture(&texture_creator).unwrap();
        canvas.clear();
        for event in event_pump.poll_iter() {
            match event {
                Event::Quit {..} |
                Event::KeyDown { keycode: Some(Keycode::Escape), .. } => {//@TODO: Clean up of threads
                    break 'running
                },
                Event::KeyDown { keycode: Some(Keycode::Space), .. } => {
                    show_samples = show_samples ^ true;
                },
                _ => {}
            }
        }
        canvas.copy(&texture,None,None).unwrap();
        canvas.present();
        ::std::thread::sleep(Duration::new(0, 1_000_000_000u32 / 2));
    }
}
