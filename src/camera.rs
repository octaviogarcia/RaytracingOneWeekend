use crate::vec3::*;
use crate::ray::*;
use crate::utils::degrees_to_radians;

#[derive(Copy,Clone)]
pub struct Camera{
    //These are all related and shouldn't really be changed individually
    //Rust sadly doesnt have const fields
    origin:            Point3,
    horizontal:        Vec3,
    vertical:          Vec3,
    lower_left_corner: Point3,
    u_of_plane: Vec3,//"Up" unit plane vector
    v_of_plane: Vec3,//"Right" unit plane vector
    w_of_plane: Vec3,//Normal unit plane vector
    lens_radius: f32,//How far is the focal point from the plane
}

impl Camera{
    pub fn world_camera(vfov_in_degrees:f32,aspect_ratio: f32) -> Self{
        return Self::new(Point3::ZERO,Point3::new(0.,0.,-1.),Vec3::new(0.,1.,0.),vfov_in_degrees,aspect_ratio,0.,1.);
    }
    //vup defines the rotation of the camera, so where your hair is (or baldspot)
    //aperture = 0 focus_dist = 1 for "Normal camera"
    pub fn new(lookfrom: Point3,lookat: Point3,vup: Vec3,vfov_in_degrees: f32,aspect_ratio: f32,aperture: f32,focus_dist: f32) -> Self{
        let vfov = degrees_to_radians(vfov_in_degrees);
        //if vfov is the angle between the bottom of the cone to the top of the cone, vfov/2 is angle the right triangle
        //if we set our plane at z=-1 then tan(vfov/2) = opposite/adjacent = opposite/1
        //"opposite" is half our height
        let height = (vfov/2.0).tan()*focus_dist;
        let viewport_height = 2.0*height;
        let viewport_width = viewport_height * aspect_ratio;

        let w_of_plane = (lookfrom - lookat).unit();
        let u_of_plane = vup.cross(w_of_plane).unit();
        let v_of_plane = w_of_plane.cross(u_of_plane);

        let o   = lookfrom;
        let h   = viewport_width  * u_of_plane;
        let v   = viewport_height * v_of_plane;

        let llc = o - h/2. - v/2. - focus_dist*w_of_plane;

        return Camera{ origin: o, horizontal: h, vertical: v, lower_left_corner: llc,
            u_of_plane: u_of_plane,v_of_plane: v_of_plane,w_of_plane: w_of_plane, lens_radius: aperture / 2.0};
    }
    pub fn get_ray(&self,u: f32,v: f32) -> Ray{
        let rand_in_lens = self.lens_radius * Vec3::rand_in_unit_disc();
        let offset = self.u_of_plane*rand_in_lens.x() + self.v_of_plane*rand_in_lens.y();
        let direction = self.lower_left_corner + u*self.horizontal + v*self.vertical - self.origin;
        return Ray::new(&(self.origin+offset),&(direction-offset).unit());
    }
}