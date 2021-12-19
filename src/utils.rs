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
pub fn normalize_color(color: &Color) -> Color{
    let r = clamp(color.x().sqrt(),0.,0.999);
    let g = clamp(color.y().sqrt(),0.,0.999);
    let b = clamp(color.z().sqrt(),0.,0.999);
    return Color::new(r,g,b);
}

/*
pub fn denormalize_color(color: &Color) -> Color{
    let r = color.x()*color.x();
    let g = color.y()*color.y();
    let b = color.z()*color.z();
    return Color::new(r,g,b);
}
#[allow(dead_code)]
pub fn write_ppm(colors: &Vec<Color>,samples_per_pixel: u32,image_width: u32,image_height: u32){
    let mut colors_str: String = "".to_owned();
    for c in colors{
        let cn = normalize_color(*c,samples_per_pixel);
        let c_str: String = format!("{} {} {}\n",256.*cn.x(),256.*cn.y(),256.*cn.z()).to_owned();
        colors_str.push_str(&c_str);
    }
    println!("P3\n{} {}\n255",image_width,image_height);
    print!("{}",colors_str);
}*/

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

#[inline]
pub fn hash1_u64(obj_id: u64) -> u64 {
    return (obj_id*456894789+348764781)%17287318477382145149;
}
#[inline]
pub fn hash2_u64(obj_id: u64) -> u64 {
    return (obj_id*56456+2345)%10520185020478678957;
}
#[inline]
pub fn hash3_u64(obj_id: u64) -> u64 {
    return (obj_id*12337+7878)%6100366985798845493;
}
#[inline]
pub fn hash4_u64(obj_id: u64) -> u64 {
    return (obj_id*7438554325+2554)%2581451885731034521;
}
#[inline]
pub fn hash5_u64(obj_id: u64) -> u64 {
    return (obj_id*12345+123123044)%2015400956511055807;
}
#[inline]
pub fn hash6_u64(obj_id: u64) -> u64 {
    return (obj_id*6373412378+12452)%8800267423223100703;
}
#[inline]
pub fn hash7_u64(obj_id: u64) -> u64 {
    return (obj_id*3453453+7874856378)%7039701875810786467;
}
#[inline]
pub fn hash8_u64(obj_id: u64) -> u64 {
    return (obj_id*999465+143)%3008457310659543551;
}
#[inline]
pub fn hash9_u64(obj_id: u64) -> u64 {
    return (obj_id*14444+111345)%5935720376112203207;
}

#[inline]//Xorshift scramble
pub fn scramble(obj_id: u64) -> u64 {//maps 0 to 0
    let mut id = obj_id;
    id ^= id << 13;
    id ^= id >>  7;
    id ^= id << 17;
    return id;
}

#[inline]
pub fn u64_to_color(id: u64) -> (u8,u8,u8) {
    let b1: u8 = ((id >>  0) & 0b11111111) as u8;
    let b2: u8 = ((id >>  8) & 0b11111111) as u8;
    let b3: u8 = ((id >> 16) & 0b11111111) as u8;
    let b4: u8 = ((id >> 24) & 0b11111111) as u8;
    let b5: u8 = ((id >> 32) & 0b11111111) as u8;
    let b6: u8 = ((id >> 40) & 0b11111111) as u8;
    let b7: u8 = ((id >> 48) & 0b11111111) as u8;
    let b8: u8 = ((id >> 56) & 0b11111111) as u8;
    return (b1 ^ b8 ^ b4,b2 ^ b5 ^ b6,b3 ^ b7);
}

#[derive(Copy,Clone)]
pub struct BloomFilter{
    pub state: u64,
}

impl BloomFilter {
    pub fn new() -> Self{ Self{state: 0} }
    pub fn set(&mut self,id: u64){
        self.state |= Self::get_hash(id);
    }
    //Returns 0 if bits aren't set, 1 if bits are set but so are others, 2 if bits are set and it's the only one
    #[allow(dead_code)]
    pub fn is_set(&self,id: u64) -> u32{
        let h = Self::get_hash(id);
        let set      = (self.state & h) == h;
        let only_set = self.state == h;
        return set as u32 + only_set as u32;
    }
    #[inline]
    fn get_hash(id: u64) -> u64{
        let mut ret = 0;
        let bit = (id != 0) as u64;//If id == 0, doesnt set anything
        ret |= bit << (hash1_u64(id) % 64);
        ret |= bit << (hash2_u64(id) % 64);
        ret |= bit << (hash3_u64(id) % 64);
        ret |= bit << (hash4_u64(id) % 64);
        ret |= bit << (hash5_u64(id) % 64);
        ret |= bit << (hash6_u64(id) % 64);
        ret |= bit << (hash7_u64(id) % 64);
        ret |= bit << (hash8_u64(id) % 64);
        ret |= bit << (hash9_u64(id) % 64);
        return ret;
    }
}