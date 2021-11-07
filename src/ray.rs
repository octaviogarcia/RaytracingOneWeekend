use crate::math::vec3::{Vec3,UnitVec3,Point3};
use crate::math::mat4x4::{Mat4x4};

#[derive(Copy,Clone)]
pub struct Ray{
    pub orig: Point3,
    pub dir: Vec3
}

impl Ray{
    pub fn new(o: &Point3,d: &UnitVec3) -> Self{
        return Ray{ orig: *o, dir: d.unit() };
    } 
    pub fn at(&self,t: f32) -> Point3{
        return self.orig + t*self.dir;
    }
    pub fn transform(&self,m: &Mat4x4) -> Ray{
        return Ray{orig: m.dot_p3(&self.orig), dir: m.dot_v3(&self.dir)};
    }
}