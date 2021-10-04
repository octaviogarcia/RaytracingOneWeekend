use crate::vec3;
use vec3::Vec3;
use vec3::Point3;

#[derive(Copy,Clone)]
pub struct Ray{
    pub orig: Point3,
    pub dir: Vec3
}

impl Ray{
    pub fn new(o: Point3,d: Vec3) -> Self{
        return Ray{ orig: o, dir: d };
    } 
    pub fn at(&self,t: f64) -> Point3{
        return self.orig + t*self.dir;
    }
}