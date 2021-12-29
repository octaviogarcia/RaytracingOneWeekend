
use crate::math::vec3::{Vec3,UnitVec3,Point3};
use crate::ray::Ray;
use crate::materials::Material;
use crate::math::mat4x4::Mat4x4;

pub trait Marched {
    fn sdf(&self,p: &Point3) -> f32;
    fn get_outward_normal(&self,p: &Point3) -> Vec3;
    fn material(&self) -> &Material;
    fn center(&self) -> &Point3;
    fn get_id(&self) -> u64;
    fn build_bounding_box(&self,look_at_inv: &Mat4x4) -> ();
    fn hit_bounding_box(&self,r: &Ray,t_min: f32,t_max: f32) -> bool;
}

pub fn get_outward_numeric_normal<M: Marched>(marched: &M,p: &Point3) -> UnitVec3{
    let eps = 0.0000001;
    let ex = Point3::new(eps, 0., 0.);
    let ey = Point3::new( 0.,eps, 0.);
    let ez = Point3::new( 0., 0.,eps);
    let x = marched.sdf(&(*p+ex)) - marched.sdf(&(*p-ex));
    let y = marched.sdf(&(*p+ey)) - marched.sdf(&(*p-ey));
    let z = marched.sdf(&(*p+ez)) - marched.sdf(&(*p-ez));
    let normal = Vec3::new(x,y,z).unit();

    //Flip the sign so always the SDF grows in the direction of the normal
    let test_ray = Ray::new(marched.center(),&normal);
    let start     = test_ray.at(0.);
    let start_val = marched.sdf(&start);
    let end     = test_ray.at(1.);
    let end_val = marched.sdf(&end);
    let sign = [-1.,1.][(end_val > start_val) as usize];//If it grows, keep the sign. Else flip it
    return normal*sign;
}

#[derive(Copy, Clone)]
pub struct MarchedSphere {
    pub center: Point3,
    pub radius: f32,
    pub material: Material
}

impl Marched for MarchedSphere {
    fn sdf(&self,p: &Point3) -> f32 {
        return (*p - self.center).length() - self.radius;
    }
    fn get_outward_normal(&self,p: &Point3) -> Vec3 {
        let normal = (*p - self.center).unit();
        return normal;
    }
    fn material(&self) -> &Material{
        return &self.material;
    }
    fn center(&self) -> &Point3{
        return &self.center;
    }
    fn get_id(&self) -> u64 { (self as *const Self) as u64 }
    fn build_bounding_box(&self,_look_at_inv: &Mat4x4) -> () {}
    fn hit_bounding_box(&self,_r: &Ray,_t_min: f32,_t_max: f32) -> bool { true }
}

#[derive(Copy, Clone)]
pub struct MarchedBox {
    pub center: Point3,
    pub sizes: Vec3,
    pub material: Material
}

impl Marched for MarchedBox {
    fn sdf(&self,p: &Point3) -> f32 {
        let q = (&(*p-self.center)).abs() - self.sizes;
        return q.max(&Vec3::ZERO).length() + q.x().max(q.y().max(q.z())).min(0.);
    }
    fn get_outward_normal(&self,p: &Point3) -> Vec3 {
        return get_outward_numeric_normal(self,p);
    }
    fn material(&self) -> &Material{
        return &self.material;
    }
    fn center(&self) -> &Point3{
        return &self.center;
    }
    fn get_id(&self) -> u64 { (self as *const Self) as u64 }
    fn build_bounding_box(&self,_look_at_inv: &Mat4x4) -> () {}
    fn hit_bounding_box(&self,_r: &Ray,_t_min: f32,_t_max: f32) -> bool { true }
}

#[derive(Copy, Clone)]
pub struct MarchedTorus {
    pub center: Point3,
    pub sizes: Vec3,//Vec2... actualy
    pub material: Material
}

impl Marched for MarchedTorus {
    fn sdf(&self,p: &Point3) -> f32 {
        let p2 = *p - self.center;
        let q = Vec3::new(Vec3::new(p2.x(),p2.z(),0.).length()-self.sizes.x(),p2.y(),0.);
        return q.length() - self.sizes.y();
    }
    fn get_outward_normal(&self,p: &Point3) -> Vec3 {
        return get_outward_numeric_normal(self,p);
    }
    fn material(&self) -> &Material{
        return &self.material;
    }
    fn center(&self) -> &Point3{
        return &self.center;
    }
    fn get_id(&self) -> u64 { (self as *const Self) as u64 }
    fn build_bounding_box(&self,_look_at_inv: &Mat4x4) -> () {}
    fn hit_bounding_box(&self,_r: &Ray,_t_min: f32,_t_max: f32) -> bool { true }
}