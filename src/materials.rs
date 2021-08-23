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

//LAMBERTIAN
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

//METAL
pub struct Metal {
    pub albedo: Color,
    pub fuzz:   f64,
}
impl Metal{
    pub fn new(a: Color) -> Self{
        return Self{albedo: a, fuzz: 0.};
    }
    pub fn new_fuzz(a: Color,f: f64) -> Self{
        return Self{albedo: a, fuzz: f};
    }
}

fn reflect(v: &Vec3,n: &Vec3) -> Vec3{
    return *v - 2.*v.dot(*n)*(*n);
}

impl Material for Metal {
    fn scatter(&self,r_in: &Ray,hr: &HitRecord) -> Option<MaterialScatterResult>{
        let reflected: Vec3 = reflect(&r_in.dir.unit(), &hr.normal);
        let new_ray = Ray::new(hr.point, reflected + self.fuzz*Vec3::rand_in_unit_sphere());
        return Some(MaterialScatterResult{attenuation: self.albedo,ray: new_ray});
    }
}