use crate::vec3;
use vec3::Vec3;
use vec3::Point3;
use vec3::Color;

use crate::ray::Ray;
use crate::materials::Material;
use std::rc::Rc;

pub struct HitRecord {
    pub point: Point3,
    pub normal: Vec3,
    pub material: Rc<dyn Material>,
    pub t: f64,
    pub front_face: bool,
}

impl HitRecord{
    pub fn set_face_normal(self: &mut Self,r: &Ray,outward_normal: &Vec3) -> () {
        self.front_face = r.dir.dot(*outward_normal) < 0.;
        if self.front_face {
            self.normal = *outward_normal;
        }
        else{
            self.normal = -*outward_normal;
        }
    }
}

pub struct Sphere {
    pub center: Point3,
    pub radius: f64,
    pub material: Rc<dyn Material>
}

pub trait Hittable {
    fn hit(&self,r: &Ray,t_min: f64,t_max: f64) -> Option<HitRecord>;
}

impl Hittable for Sphere {
    fn hit(&self,r: &Ray,t_min: f64,t_max: f64) -> Option<HitRecord> {
        let oc = r.orig - self.center;
        let a = r.dir.length_squared();
        let half_b = oc.dot(r.dir);
        let c = oc.length_squared() - self.radius*self.radius;
        let discriminant = half_b*half_b - a*c;
        if discriminant < 0. {
            return None;
        }
        let sqrtd = discriminant.sqrt();
        let mut root = (-half_b - sqrtd)/a;
        if root < t_min || root > t_max{
            root = (-half_b + sqrtd)/a;
            if root < t_min || root > t_max{
                return None;
            }
        }
        let point = r.at(root);
        let outward_normal = (point - self.center)/self.radius;
        let mut rec = HitRecord{t: root,point: point,normal: Vec3::ZERO, front_face: false,material: self.material.clone()};
        rec.set_face_normal(r,&outward_normal);
        return Some(rec);
    }
}

pub struct HittableList{
    pub objects: Vec<Box<dyn Hittable>>
}

impl HittableList{
    pub fn new() -> Self {
        return HittableList{objects: Vec::new()};
    }
    pub fn add(&mut self,obj: Box<dyn Hittable>) -> () {
        self.objects.push(obj);
    }
    pub fn clear(&mut self) -> () {
        self.objects.clear();
    }
    pub fn hit(&self,r: &Ray,t_min: f64,t_max: f64) -> Option<HitRecord> {
        let mut closest_so_far = t_max;
        let mut rec: Option<HitRecord>  = None;
        for obj in &self.objects{
            match obj.hit(r,t_min,closest_so_far) {
                Some(hr) => {
                    if hr.t < closest_so_far {
                        closest_so_far = hr.t;
                        rec = Some(hr);
                    }
                }
                None => {}
            }
        }
        return rec;
    }
}