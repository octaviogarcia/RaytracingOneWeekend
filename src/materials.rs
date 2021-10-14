use crate::ray::Ray;
use crate::vec3::*;
use crate::hits::HitRecord;
use crate::utils::min;
use crate::utils::MyRandom;

pub struct MaterialScatterResult {
    pub attenuation: Color,
    pub ray: Ray,
}

#[derive(Copy, Clone, PartialEq)]
pub enum MaterialType {
    LAMBERTIAN,
    METAL,
    DIELECTRIC,
}

#[derive(Copy, Clone)]
pub struct Material {//Used in:
    pub albedo: Color,//Lambertian, Metal
    pub fuzz: f32,//Metal
    pub ior: f32,//Dielectric
    pub mat_type: MaterialType, //Tag
}

impl Material{
    pub fn new_lambertian(albedo: Color) -> Self{
        return Self{albedo: albedo,fuzz: 0.,ior: 0.,mat_type: MaterialType::LAMBERTIAN};
    }
    pub fn new_metal(albedo: Color) -> Self{
        return Self{albedo: albedo,fuzz: 0.,ior: 0.,mat_type: MaterialType::METAL};
    }
    pub fn new_metal_fuzz(albedo: Color,fuzz: f32) -> Self{
        return Self{albedo: albedo,fuzz: fuzz,ior: 0.,mat_type: MaterialType::METAL};
    }
    pub fn new_dielectric(index_of_refraction: f32) -> Self{
        return Self{albedo: Color::ZERO,fuzz: 0.,ior: index_of_refraction,mat_type: MaterialType::DIELECTRIC};
    }
    pub fn scatter(&self,r_in: &Ray,hr: &HitRecord) -> MaterialScatterResult {
        match &self.mat_type{
            MaterialType::LAMBERTIAN => {
                return self.scatter_lambertian(r_in,hr);
            }
            MaterialType::METAL => {
                return self.scatter_metal(r_in,hr);
            }
            MaterialType::DIELECTRIC => {
                return self.scatter_dielectric(r_in,hr);
            }
        }
    }
    pub fn scatter_lambertian(&self,_r_in: &Ray,hr: &HitRecord) -> MaterialScatterResult {
        //let new_dir = hr.normal + Vec3::rand_in_unit_sphere();
        let mut new_dir = hr.normal + Vec3::rand_unit_vector();
        if new_dir.near_zero() {//If by offchance we make it zero, just use the normal
            new_dir = hr.normal;
        }
        //let new_dir = Vec3::rand_in_hemisphere(&hr.normal);
        let new_ray = Ray::new(hr.point, new_dir);
        return MaterialScatterResult{attenuation: self.albedo,ray: new_ray};
    }
    pub fn scatter_metal(&self,r_in: &Ray,hr: &HitRecord) -> MaterialScatterResult {
        let reflected: Vec3 = reflect(&r_in.dir.unit(), &hr.normal);
        let new_ray = Ray::new(hr.point, reflected + self.fuzz*Vec3::rand_in_unit_sphere());
        return MaterialScatterResult{attenuation: self.albedo,ray: new_ray};
    }

    pub fn scatter_dielectric(&self,r_in: &Ray,hr: &HitRecord) -> MaterialScatterResult {
        let dir_unit = r_in.dir.unit();
        let front_face = dir_unit.dot(hr.normal) < 0.;
        let refraction_ratio: f32;
        let normal: Vec3;
        if front_face {
            refraction_ratio = 1.0/self.ior;
            normal = hr.normal;
        }
        else{
            refraction_ratio = self.ior;
            normal = -hr.normal;
        }

        let cos_theta = min((-dir_unit).dot(normal),1.0);
        let sin_theta = (1.0 - cos_theta*cos_theta).sqrt();
        let cannot_refract = (refraction_ratio * sin_theta) > 1.0;
        let reflect_by_reflectante = reflectance(cos_theta, refraction_ratio) > f32::rand();
        let new_dir: Vec3;
        if cannot_refract || reflect_by_reflectante {//Reflect
            new_dir = reflect(&dir_unit, &normal);
        }
        else {
            new_dir = refract(&dir_unit,&normal,refraction_ratio);
        }

        return MaterialScatterResult{attenuation:  Color::new(1.0,1.0,1.0),ray: Ray::new(hr.point,new_dir.unit())};
    }
}

fn reflect(v: &Vec3,n: &Vec3) -> Vec3{
    return *v - 2.*v.dot(*n)*(*n);
}

fn refract(uv: &Vec3,n: &Vec3,etai_over_etat: f32) -> Vec3 {
    let cos_theta = min((-*uv).dot(*n),1.0);
    let r_out_perp = etai_over_etat * ((*uv) + cos_theta*(*n));
    let aux = -(1.0 - r_out_perp.length_squared()).abs().sqrt();
    let r_out_parallel = aux*(*n);
    return r_out_perp + r_out_parallel;
}

fn reflectance(cos: f32, ref_idx: f32) -> f32{
    // Use Schlick's approximation for reflectance.
    let r0 = (1.-ref_idx)/(1.+ref_idx);
    let r0_2 = r0*r0;
    let cos_5 = (1.-cos)*(1.-cos)*(1.-cos)*(1.-cos)*(1.-cos);
    return r0_2 + (1.-r0_2)*cos_5;
}