use crate::ray::Ray;
use crate::vec3::*;
use crate::hits::HitRecord;
use crate::utils::min;
use crate::utils::MyRandom;

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

fn refract(uv: &Vec3,n: &Vec3,etai_over_etat: f64) -> Vec3 {
    let cos_theta = min((-*uv).dot(*n),1.0);
    let r_out_perp = etai_over_etat * ((*uv) + cos_theta*(*n));
    let aux = -(1.0 - r_out_perp.length_squared()).abs().sqrt();
    let r_out_parallel = aux*(*n);
    return r_out_perp + r_out_parallel;
}
fn reflectance(cos: f64, ref_idx: f64) -> f64{
    // Use Schlick's approximation for reflectance.
    let r0 = (1.-ref_idx)/(1.+ref_idx);
    let r0_2 = r0*r0;
    let cos_5 = (1.-cos)*(1.-cos)*(1.-cos)*(1.-cos)*(1.-cos);
    return r0_2 + (1.-r0_2)*cos_5;
}

pub struct Dieletric{
    pub ior: f64,
}
impl Dieletric{
    pub fn new(index_of_refraction: f64) -> Self{
        return Self{ior: index_of_refraction};
    }
}
impl Material for Dieletric {
    fn scatter(&self,r_in: &Ray,hr: &HitRecord) -> Option<MaterialScatterResult>{
        let refraction_ratio: f64;
        if hr.front_face {
            refraction_ratio = 1.0 / self.ior;
        }
        else{
            refraction_ratio = self.ior;
        }

        let cos_theta = min((-r_in.dir.unit()).dot(hr.normal),1.0);
        let sin_theta = (1.0 - cos_theta*cos_theta).sqrt();
        let cannot_refract = (refraction_ratio * sin_theta) > 1.0;
        let reflect_by_reflectante = reflectance(cos_theta, refraction_ratio) > f64::rand();
        let new_dir: Vec3;
        if cannot_refract || reflect_by_reflectante {//Reflect
            new_dir = reflect(&r_in.dir.unit(), &hr.normal);
        }
        else {
            new_dir = refract(&r_in.dir.unit(),&hr.normal,refraction_ratio);
        }
        let new_ray   = Ray::new(hr.point,new_dir);
        let attenuation = Color::new(1.0,1.0,1.0);
        return Some(MaterialScatterResult{attenuation: attenuation,ray: new_ray});
    }
}