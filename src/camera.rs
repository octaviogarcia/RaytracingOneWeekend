use crate::vec3::*;
use crate::ray::*;
use crate::utils::degrees_to_radians;

pub struct Camera{
    pub origin:            Point3,
    pub horizontal:        Vec3,
    pub vertical:          Vec3,
    pub lower_left_corner: Point3,
}

impl Camera{
    pub fn new(vfov_in_degrees: f64,aspect_ratio: f64) -> Self{
        let vfov = degrees_to_radians(vfov_in_degrees);
        //if vfov is the angle between the bottom of the cone to the top of the cone, vfov/2 is angle the right triangle
        //if we set our plane at z=-1 then tan(vfov/2) = opposite/adjacent = opposite/1
        //"opposite" is half our height
        let h = (vfov/2.0).tan();
        let viewport_height = 2.0*h;
        let viewport_width = viewport_height * aspect_ratio;
        let focal_length = 1.0;

        let h   = Vec3::new(viewport_width,0.,0.);
        let v   = Vec3::new(0.,viewport_height,0.);
        let fl  = Vec3::new(0.,0.,focal_length);
        let o   = Point3::ZERO;//Assumes origin at zero, +y up, +x right,-z inside the screen (RH rule)
        let llc = o - h/2. - v/2. - fl;
        
        return Camera{ origin: o, horizontal: h, vertical: v, lower_left_corner: llc };
    }
    pub fn get_ray(&self,u: f64,v: f64) -> Ray{
        let direction = self.lower_left_corner + u*self.horizontal + v*self.vertical - self.origin;
        return Ray::new(self.origin,direction);
    }
}