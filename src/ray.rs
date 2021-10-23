use crate::math::vec3::{Vec3,UnitVec3,Point3};

#[derive(Copy,Clone)]
pub struct Ray{
    pub orig: Point3,
    pub dir: UnitVec3
}

impl Ray{
    pub fn new(o: &Point3,d: &Vec3) -> Self{
        return Ray{ orig: *o, dir: d.unit() };
    } 
    pub fn at(&self,t: f32) -> Point3{
        return self.orig + t*self.dir;
    }
}