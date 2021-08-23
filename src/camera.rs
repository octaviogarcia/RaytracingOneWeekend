use crate::vec3;
use vec3::*;

use crate::ray;
use ray::*;

pub struct Camera{
    pub origin:            Point3,
    pub horizontal:        Vec3,
    pub vertical:          Vec3,
    pub lower_left_corner: Point3,
}

impl Camera{
    pub fn new(viewport_height: f64,viewport_width: f64,focal_length: f64) -> Self{
        let o   = Point3::ZERO;//Assumes origin at zero, +y up, +x right,-z inside the screen (RH rule)
        let h   = Vec3::new(viewport_width,0.,0.);
        let v   = Vec3::new(0.,viewport_height,0.);
        let fl  = Vec3::new(0.,0.,focal_length);
        let llc = o - h/2. - v/2. - fl;
        return Camera{ origin: o, horizontal: h, vertical: v, lower_left_corner: llc };
    }
    pub fn get_ray(&self,u: f64,v: f64) -> Ray{
        let direction = self.lower_left_corner + u*self.horizontal + v*self.vertical - self.origin;
        return Ray::new(self.origin,direction);
    }
}