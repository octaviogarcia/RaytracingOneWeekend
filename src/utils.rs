use crate::math::vec3;
use vec3::Color;

pub fn min(x: f32,y: f32) -> f32{
    if x<y{
        return x;
    }
    return y;
}
pub fn max(x: f32,y: f32) -> f32{
    if x<y{
        return y;
    }
    return x;
}
pub fn abs(x: f32) -> f32{
    let sign = x > 0.;
    let map = [-1.,1.];
    return x*map[sign as usize];
}

pub fn clamp(x: f32,fmin: f32,fmax: f32) -> f32{
    return min(max(x,fmin),fmax);
}

//t = 0 -> c1, t = 1 -> c2
pub fn lerp(t: f32,c1: Color,c2: Color) -> Color{
    return (1.0-t)*c1 + t*c2;
}
pub fn normalize_color(color: Color,samples_per_pixel: u32) -> Color{
    let scale = 1.0/(samples_per_pixel as f32);
    let r = clamp((color.x() * scale).sqrt(),0.,0.999);
    let g = clamp((color.y() * scale).sqrt(),0.,0.999);
    let b = clamp((color.z() * scale).sqrt(),0.,0.999);
    return Color::new(r,g,b);
}
pub fn denormalize_color(color: Color,samples_per_pixel: u32) -> Color{
    let scale = samples_per_pixel as f32;
    let r = clamp(color.x()*color.x()*scale,0.,scale);
    let g = clamp(color.y()*color.y()*scale,0.,scale);
    let b = clamp(color.z()*color.z()*scale,0.,scale);
    return Color::new(r,g,b);
}

pub fn write_ppm(colors: &Vec<Color>,samples_per_pixel: u32,image_width: u32,image_height: u32){
    let mut colors_str: String = "".to_owned();
    for c in colors{
        let cn = normalize_color(*c,samples_per_pixel);
        let c_str: String = format!("{} {} {}\n",256.*cn.x(),256.*cn.y(),256.*cn.z()).to_owned();
        colors_str.push_str(&c_str);
    }
    println!("P3\n{} {}\n255",image_width,image_height);
    print!("{}",colors_str);
}

use rand::Rng;
pub trait MyRandom{
    fn rand() -> Self;
    fn rand_range(fmin: f32,fmax: f32) -> Self;
}
impl MyRandom for f32{
    fn rand() -> Self{ rand::thread_rng().gen() }
    fn rand_range(min: f32,max: f32) -> f32{ Self::rand()*(max-min) + min }
}

pub const PI:  f32 = 3.1415926535897932385;
pub const INF: f32 = f32::INFINITY;

pub fn degrees_to_radians(degrees: f32) -> f32{
    return degrees * PI / 180.0;
}

pub fn round_n(f: f64,n: u32) -> f64 {
    let mut mul = 1.0;
    for _n in 0..n{
        mul *= 10.0;
    }
    return (f*mul).round()/mul;
}