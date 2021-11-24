use rand::prelude::*;
use crate::math::vec3::*;
use crate::hits::*;
use crate::camera::*;
use crate::ray::*;
use std::sync::atomic::{AtomicU64, Ordering};
use crate::utils::{lerp,MyRandom};

#[derive(Copy,Clone)]
pub struct ColorsBox {
    pub colors: *mut Vec<Color>,
    pub samples: *mut Vec<u32>,
    pub true_samples: *mut Vec<u32>,
}
unsafe impl Send for ColorsBox{}

struct Stats{//https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm
    n: f32,
    sum: Color,
    m2: Color,
    avg: Color
}

impl Stats{
    pub fn new() -> Self {
        Self{sum:Color::ZERO,m2:Color::ZERO,n:0.,avg:Color::ZERO}
    }
    #[inline]
    pub fn add(&mut self,x: &Color) -> f32{
        let old_avg = self.avg;
        self.sum   += x;
        self.n     += 1.;
        self.avg    = self.sum / self.n;
        self.m2   +=  (*x-old_avg)*(*x-self.avg);
        let var    = self.m2 / (self.n-1.);
        return ((*x-self.avg)/var.sqrt()).abs().max_val();
    }
}

struct ThreadPixels{
    //Assigned pixels to the threads
    pub indexes: Vec<usize>,
    pub len: usize,
    //When rendering, here are added pixels that need to be rendered next pass
    //After a pass is done, one is expected to swap() these 2 buffers
    pub backbuff_indexes: Vec<usize>,
    pub backbuff_len: usize,
    //Useless run per pixel, past a threshold we stop raycasting it
    pub useless_runs: Vec<u32>,
}

impl ThreadPixels{
    pub fn new(expected_max_size: usize) -> Self{
        Self{indexes: Vec::with_capacity(expected_max_size),len: 0,
             backbuff_indexes: Vec::with_capacity(expected_max_size),backbuff_len: 0,
             useless_runs: Vec::with_capacity(expected_max_size)}
    }
    pub fn push(&mut self,idx: usize){
        self.indexes.push(idx);
        self.len += 1;
        self.backbuff_indexes.push(999999999);//Add garbage just to expand if len overflows expected_max_size
        self.useless_runs.push(0);//Init to 0
    }
    pub fn shrink_to_fit(&mut self){
        self.indexes.shrink_to_fit();
        self.backbuff_indexes.shrink_to_fit();
        self.useless_runs.shrink_to_fit();
    }
    pub fn swap_buffers(&mut self){
        ::std::mem::swap(&mut self.indexes,&mut self.backbuff_indexes);
        self.len = self.backbuff_len;
        self.backbuff_len = 0;
    }
    pub fn add_useless_run(&mut self,is_useless: bool,i: usize) -> bool{//True if useless runs go past the threshold
        self.useless_runs[i] += is_useless as u32;
        self.useless_runs[i] *= is_useless as u32;
        let mur = self.useless_runs[i] == 10;//@TODO: make more configurable
        //If the pixel render is "not useless", keep adding to the back buffer to draw in the next iteration
        self.backbuff_indexes[self.backbuff_len] = self.indexes[i];
        self.backbuff_len+=(!mur) as usize;
        return mur;
    }
}

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

pub fn render(camera: &Camera,world: &FrozenHittableList,max_depth: u32,tmin: f32,tmax: f32,
    samples_per_pixel: u32,image_width: u32,image_height: u32,
    colors_box: ColorsBox,
    tid: u32,assigned_thread: &Vec<u32>,samples_atom: &AtomicU64)
{
    let image_width_f  = image_width as f32;
    let image_height_f = image_height as f32;
    let image_size = (image_width*image_height) as usize;

    let mut thread_pixels = ThreadPixels::new(image_size);
    let mut stats: Vec<Stats> = Vec::with_capacity(image_size);
    for pos in 0..image_size {
        if assigned_thread[pos] == tid {
            thread_pixels.push(pos);
            stats.push(Stats::new());
        }
    }
    thread_pixels.shrink_to_fit();
    stats.shrink_to_fit();

    //Construct a blue noise wannabe with a low discrepancy random number
    //https://en.wikipedia.org/wiki/Low-discrepancy_sequence#Construction_of_low-discrepancy_sequences
    let jitters: Vec<(f32,f32)> = {
        let mut ret: Vec<(f32,f32)> = Vec::with_capacity(samples_per_pixel as usize);
        for s in 0..samples_per_pixel{
            let sdiv = s / 2;
            let jitteri = (sdiv&1) as f32;//mod 2
            let jitterj = (s&1) as f32;//mod 2
            ret.push((jitteri,jitterj));
        }//produces (0,0),(0,1),(1,0),(1,1),(0,0),...
        ret.shuffle(&mut rand::thread_rng());
        ret
    };

    //If samples_per_pixel > 30, in the first 30 runs we simply render and gain statistics
    //after that, we do the optimization
    let minimal_stats_runs = samples_per_pixel.min(30); 

    for _sample in 0..samples_per_pixel{
        //posidx is a double indirection... indexes[pos_idx] is the pixel index in the main memory buffer
        for pos_idx in 0..thread_pixels.len{
            let idx = thread_pixels.indexes[pos_idx];
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
                (*(colors_box.colors))[idx]       += pixel_color; 
                (*(colors_box.samples))[idx]      += 1;
                (*(colors_box.true_samples))[idx] += 1;

                let curr_color = (*(colors_box.colors))[idx]/(curr_samples as f32 + 1.);
                let stat = &mut stats[pos_idx];
                let max_abs_z = stat.add(&curr_color);
                //We set the threshold at 2 stddevs
                //If we are in the initial runs, always "adds" false, disabling per se the mechanism
                let mur     = thread_pixels.add_useless_run(max_abs_z < 1.5 && curr_samples > minimal_stats_runs,pos_idx) as u32;
                let mur_neg = 1 - mur;
                //When enough useless runs go by we simply set it to the max (samples_per_pixel), else increment
                let aux_samples = mur*samples_per_pixel+mur_neg*(curr_samples+1);
                (*(colors_box.samples))[idx] = aux_samples;
                (*(colors_box.colors))[idx]  = (aux_samples as f32)*curr_color;

                //Inform sample is done to Log Thread
                //We fetch add left over samples, or 1
                samples_atom.fetch_add((aux_samples-curr_samples) as u64,Ordering::Relaxed);
            }
        }
        thread_pixels.swap_buffers();
    }
}