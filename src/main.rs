extern crate num_cpus;
extern crate sdl2;

mod math;
use math::vec3::*;

mod utils;
use utils::*;

mod ray;
use ray::*;

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

use sdl2::event::Event;
use sdl2::keyboard::Keycode;
use std::time::Duration;
use rand::prelude::*;

#[derive(Copy,Clone)]
struct ColorsBox {
    pub colors: *mut Vec<Color>,
    pub samples: *mut Vec<u32>,
    pub true_samples: *mut Vec<u32>,
}
unsafe impl Send for ColorsBox{}

fn ray_color(r: &Ray,world: &FrozenHittableList, depth: u32,tmin: f32,tmax: f32) -> Color{
    let mut curr_color = Color::new(1.,1.,1.);
    let mut curr_ray: Ray = *r;
    for _i in 0..depth{
        match world.hit(&curr_ray,tmin,tmax) {
            Some(hr) => {
                let rslt = hr.material.scatter(r,&hr);
                curr_color *= rslt.attenuation;
                curr_ray = rslt.ray;
            },
            None => {
                let t: f32 = 0.5*(r.dir.y() + 1.0);
                let lerped_sky_color = lerp(t,Color::new(1.0,1.0,1.0),Color::new(0.5,0.7,1.0));
                return curr_color*lerped_sky_color;
            }
        }
    }
    return -Color::ZERO;//If we run out of depth return -black
}


fn random_scene() -> HittableList{
    let mut world = HittableList::new();
    let mat_ground = Material::new_lambertian(Color::new(0.5,0.5,0.5));
    //world+=&MarchedSphere{center: Point3::new(0., -1000.,0.), radius: 1000.0, material: mat_ground};
    world+=&Sphere{center: Point3::new(0., -1000.,0.), radius: 1000.0, material: mat_ground};
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
                world+=&Sphere{center: center, radius: 0.2, material: sphere_mat};
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
        world+=&Cube{center: Point3::new(4.,1.,0.), radius: 0.5, material: mat};
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

fn draw(camera: &Camera,world: &FrozenHittableList,max_depth: u32,tmin: f32,tmax: f32,
    samples_per_pixel: u32,image_width: u32,image_height: u32,
    colors_box: ColorsBox,
    tid: u32,assigned_thread: &Vec<u32>,samples_atom: &AtomicU64)
{
    let image_width_f  = image_width as f32;
    let image_height_f = image_height as f32;
    let image_size = (image_width*image_height) as usize;

    let mut thread_pixels:       Vec<usize> = Vec::with_capacity(image_size);
    let mut thread_pixels_len:  usize;
    let mut thread_pixels_back:  Vec<usize> = Vec::with_capacity(image_size);
    let mut thread_pixels_back_len:  usize = 0;
    let mut thread_pixels_useless_runs: Vec<u32>   = Vec::with_capacity(image_size);
    for pos in 0..image_size {
        if assigned_thread[pos] == tid {
            thread_pixels.push(pos);
            thread_pixels_back.push(9999999);
            thread_pixels_useless_runs.push(0);
        }
    }

    thread_pixels.shrink_to_fit();
    thread_pixels_len = thread_pixels.len();
    thread_pixels_back.shrink_to_fit();
    thread_pixels_useless_runs.shrink_to_fit();

    //Flip it by thread odness so we don't give priority to starting pixels while drawing
    if (tid % 2) == 1 { thread_pixels = thread_pixels.into_iter().rev().collect(); }

    const EPS: f32 = 0.01; 
    const MAX_USELESS_RUNS: u32 = 5;


    //Construct a blue noise wannabe with a low discrepancy random number
    //https://en.wikipedia.org/wiki/Low-discrepancy_sequence#Construction_of_low-discrepancy_sequences
    let jitters: Vec<(f32,f32)> = {
        let mut ret: Vec<(f32,f32)> = Vec::with_capacity(samples_per_pixel as usize);
        for s in 0..samples_per_pixel{
            let div: u32 = s / 2;
            let mmod: u32 = s - div*2;//@TODO: Pregenerate jitters and shuffle por proper
            let jitteri = (div&1) as f32;//mod 2
            let jitterj = mmod as f32;
            ret.push((jitteri,jitterj));
        }
        ret.shuffle(&mut rand::thread_rng());
        ret
    };

    for _sample in 0..samples_per_pixel{
        for pos_idx in 0..thread_pixels_len{
            let idx = thread_pixels[pos_idx];
            let curr_samples = unsafe { (*colors_box.samples)[idx] };
            //Should never happen since we upkeep undone pixels with a backbuffer
            //assert!(curr_samples < samples_per_pixel);

            let line = (idx as u32) / image_width;
            let col  = (idx as u32) - image_width*line;
            let j_f = line as f32;
            let i_f = col as f32;

            let i_rand = (f32::rand() + jitters[curr_samples as usize].0)/2.;
            let j_rand = (f32::rand() + jitters[curr_samples as usize].1)/2.;
            let u = (i_f+i_rand)/(image_width_f-1.);
            let v = (j_f+j_rand)/(image_height_f-1.);
            let ray = camera.get_ray(u,1.0-v);
            let pixel_color = ray_color(&ray,&world,max_depth,tmin,tmax);

            unsafe {
                let old_color = (*(colors_box.colors))[idx]/(curr_samples as f32);

                (*(colors_box.colors))[idx]       += pixel_color; 
                (*(colors_box.samples))[idx]      += 1;
                (*(colors_box.true_samples))[idx] += 1;

                let curr_color = (*(colors_box.colors))[idx]/(curr_samples as f32 + 1.);
                let var = ((curr_color/old_color) - Color::new(1.,1.,1.)).abs();
                let stats = max(max(var.x(),var.y()),var.z());
                
                thread_pixels_useless_runs[pos_idx] += (stats < EPS) as u32;
                thread_pixels_useless_runs[pos_idx] *= (stats < EPS) as u32;

                let mur     =  (thread_pixels_useless_runs[pos_idx] == MAX_USELESS_RUNS) as u32;
                let mur_neg = (!(thread_pixels_useless_runs[pos_idx] == MAX_USELESS_RUNS)) as u32;
                let aux_samples = mur*samples_per_pixel+mur_neg*(curr_samples+1);
                (*(colors_box.samples))[idx] = aux_samples;
                (*(colors_box.colors))[idx]  = (aux_samples as f32)*curr_color;
                samples_atom.fetch_add((aux_samples-curr_samples-1) as u64,Ordering::Relaxed);//-1 cause we fetch_add downthere

                //If the pixel render is "not useless", keep adding to the back buffer to draw in the next iteration
                thread_pixels_back[thread_pixels_back_len] = idx;
                thread_pixels_back_len+=mur_neg as usize;
            }

            //Inform sample is done to Log Thread
            samples_atom.fetch_add(1,Ordering::Relaxed);
        }
        ::std::mem::swap(&mut thread_pixels,&mut thread_pixels_back);
        thread_pixels_len = thread_pixels_back_len;
        thread_pixels_back_len = 0;
    }
}

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
            loop {
                let progress = smpls_atom.load(Ordering::Relaxed);
                print_progress((progress as f64)/total_samples_f);
                if total_samples == progress{ 
                    print_progress(1.0);
                    return;
                }
                ::std::thread::sleep(Duration::new(0, 1_000_000_000u32 / 2));
            }
        });
    }

    let num_threads = num_cpus::get() as u32 - 1;

    let mut assigned_thread: Vec<u32> = Vec::with_capacity(image_size as usize);
    for cidx in 0..image_size{
        assigned_thread.push(cidx%num_threads);
    }
    assigned_thread.shuffle(&mut rand::thread_rng());
    assigned_thread.shrink_to_fit();
    let arc_assigned_thread = Arc::new(assigned_thread);

    let colors_box = ColorsBox{colors: &mut vec!(Color::ZERO;image_size as usize),
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
            return draw(&cam,&wrld.freeze(&cam),max_depth,tmin,tmax,
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
    let mut texture = texture_creator
    .create_texture_target(texture_creator.default_pixel_format(), image_width as u32, image_height as u32)
    .unwrap();

    let mut show_samples = false;

    'running: loop {
        canvas.with_texture_canvas(&mut texture, |texture_canvas| {
            if show_samples {
                let mut max_samples = 1.0;
                for pos in 0..(image_width*image_height){
                    if (true_samples[pos as usize] as f32) > max_samples {
                        max_samples = true_samples[pos as usize] as f32;
                    }
                }
                for pos in 0..(image_width*image_height){
                    let smpls = (true_samples[pos as usize] as f32)/max_samples;
                    let c = normalize_color(Color::new(smpls,smpls,smpls),1);
                    texture_canvas.set_draw_color(sdl2::pixels::Color::RGB((c.x()*256.0) as u8,(c.y()*256.0) as u8,(c.z()*256.0) as u8));
                    let y = pos / image_width;
                    let x = pos - y*image_width;
                    texture_canvas.draw_point(sdl2::rect::Point::new(x as i32, y as i32)).unwrap();
                }
                return;
            }

            for pos in 0..(image_width*image_height){
                let smpls = samples[pos as usize];
                let c = normalize_color(colors[pos as usize],smpls);
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
