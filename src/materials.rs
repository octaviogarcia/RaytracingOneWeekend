use crate::ray::Ray;
use crate::vec3::*;
use crate::hits::HitRecord;

pub struct MaterialScatterResult {
    pub attenuation: Color,
    pub ray: Ray,
}

pub trait Material {
    fn scatter(&self,r_in: &Ray,rec: &HitRecord) -> Option<MaterialScatterResult>;
}

pub struct Lambertian {
    pub albedo: Color,
}

impl Lambertian{
    pub fn new(a: Color) -> Self{
        return Self{albedo: a};
    }
}
impl Material for Lambertian {
    fn scatter(&self,r_in: &Ray,hr: &HitRecord) -> Option<MaterialScatterResult>{
        //let new_dir = hr.normal + Vec3::rand_in_unit_sphere();
        let mut new_dir = hr.normal + Vec3::rand_unit_vector();
        if new_dir.near_zero() {//If by offchance we make it zero, just use the normal
            new_dir = hr.normal;
        }
        //let new_dir = Vec3::rand_in_hemisphere(&hr.normal);
        let new_ray = Ray::new(hr.point, new_dir);
        return Some(MaterialScatterResult{attenuation: self.albedo,ray: new_ray});
    }
}