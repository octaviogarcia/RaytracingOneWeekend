use crate::vec3;
use vec3::Color;

pub fn min(x: f64,y: f64) -> f64{
    if x<y{
        return x;
    }
    return y;
}
pub fn max(x: f64,y: f64) -> f64{
    if x<y{
        return y;
    }
    return x;
}

pub fn clamp(x: f64,fmin: f64,fmax: f64) -> f64{
    return min(max(x,fmin),fmax);
}

//t = 0 -> c1, t = 1 -> c2
pub fn lerp(t: f64,c1: Color,c2: Color) -> Color{
    return (1.0-t)*c1 + t*c2;
}
pub fn write_color(&color: &Color,samples_per_pixel: u64){
    let scale = 1.0/(samples_per_pixel as f64);
    let r = clamp((color.x() * scale).sqrt(),0.,0.999);
    let g = clamp((color.y() * scale).sqrt(),0.,0.999);
    let b = clamp((color.z() * scale).sqrt(),0.,0.999);
    println!("{} {} {}",256.*r,256.*g,256.*b);
}

use rand::Rng;
pub trait MyRandom{
    fn rand() -> Self;
    fn rand_range(fmin: f64,fmax: f64) -> Self;
}
impl MyRandom for f64{
    fn rand() -> Self{ rand::thread_rng().gen() }
    fn rand_range(min: f64,max: f64) -> f64{ Self::rand()*(max-min) + min }
}

pub const PI:  f64 = 3.1415926535897932385;
pub const INF: f64 = f64::INFINITY;

pub fn degrees_to_radians(degrees: f64) -> f64{
    return degrees * PI / 180.0;
}
