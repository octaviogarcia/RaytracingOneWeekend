use crate::vec3;
use vec3::Vec3;
use vec3::Point3;
use crate::utils::INF;

use crate::ray::Ray;
use crate::materials::Material;

pub struct HitRecord {
    pub point: Point3,
    pub normal: Vec3,
    pub material: Material,
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

#[derive(Copy, Clone)]
pub struct Sphere {
    pub center: Point3,
    pub radius: f64,
    pub material: Material
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
        //Maybe its faster to send some sort of reference/pointer to material? Probably not, since its so small
        let mut rec = HitRecord{t: root,point: point,normal: Vec3::ZERO, front_face: false,material: self.material};
        rec.set_face_normal(r,&outward_normal);
        return Some(rec);
    }
}

#[derive(Copy, Clone)]
pub struct MarchedSphere {
    pub center: Point3,
    pub radius: f64,
    pub material: Material
}

pub trait Marched {
    fn sdf(&self,p: &Point3) -> f64;
    fn get_outward_normal(&self,p: &Point3) -> Vec3;
    fn material(&self) -> &Material;
}

impl Marched for MarchedSphere {
    fn sdf(&self,p: &Point3) -> f64 {
        return (*p - self.center).length() - self.radius;
    }
    fn get_outward_normal(&self,p: &Point3) -> Vec3 {
        return (*p - self.center).unit();
    }
    fn material(&self) -> &Material{
        return &self.material;
    }
}

pub struct HittableList{
    pub spheres: Vec<Sphere>,
    pub objects: Vec<Box<dyn Hittable + Send + Sync>>,
    pub marched_spheres: Vec<MarchedSphere>,
    pub marched_objects: Vec<Box<dyn Marched + Send + Sync>>,
}

impl HittableList{
    pub fn new() -> Self {
        return HittableList{spheres: Vec::new(),marched_spheres: Vec::new(),objects: Vec::new(),marched_objects: Vec::new()};
    }
    pub fn add_sphere(&mut self,obj: &Sphere) -> () {
        self.spheres.push(*obj);
    }
    pub fn add_marched_sphere(&mut self,obj: &MarchedSphere) -> () {
        self.marched_spheres.push(*obj);
    }
    pub fn add(&mut self,obj: Box<dyn Hittable + Send + Sync>) -> () {
        self.objects.push(obj);
    }
    pub fn add_marched(&mut self,obj: Box<dyn Marched + Send + Sync>) -> () {
        self.marched_objects.push(obj);
    }
    pub fn clear(&mut self) -> () {
        self.objects.clear();
    }

    pub fn hit(&self,r: &Ray,t_min: f64,t_max: f64) -> Option<HitRecord> {
        //Ray tracing section
        let mut closest_so_far = t_max;
        let mut rec: Option<HitRecord>  = None;
        for obj in &self.spheres{
            if let Some(hr) = obj.hit(r,t_min,closest_so_far) {
                if hr.t < closest_so_far {
                    closest_so_far = hr.t;
                    rec = Some(hr);
                }
            }
        }
        for obj in &self.objects{
            if let Some(hr) = obj.hit(r,t_min,closest_so_far) {
                if hr.t < closest_so_far {
                    closest_so_far = hr.t;
                    rec = Some(hr);
                }
            }
        }
        //Ray marching section
        const HIT_SIZE: f64 = 0.0001;
        let mut t = t_min;   
        'raymarch: while t < t_max{
            let p = r.at(t);
            let (d,normal,material) = self.get_closest_distance_normal_material(&p);
            match material {
                None => {
                    break 'raymarch;
                }
                Some(m) => {
                    if d < HIT_SIZE {//We hit something
                        if t < closest_so_far {//Its better than what we got raytracing
                            let mut hr = HitRecord{t: t,point: p,normal: Vec3::ZERO, front_face: false,material: m};
                            hr.set_face_normal(r,&normal);
                            rec = Some(hr);
                        }
                        //Else, We hit something but its not good enough
                        break 'raymarch;
                    }
                    //Move forward
                    else { t += d; }//Don't think this is correct
                }
            }
        }
        return rec;
    }

    pub fn get_closest_distance_normal_material(&self,p: &Point3) -> (f64,Vec3,Option<Material>){
        let mut max_dis = INF;
        let mut normal = Vec3::ZERO;
        let mut material = None;
        for obj in &self.marched_spheres {
            let d = obj.sdf(p);//Not an actual vtable call, just a normal fast function call
            if d < max_dis {
                max_dis = d;
                normal = obj.get_outward_normal(p);//Not an actual vtable call, just a normal fast function call
                material = Some(obj.material);
            }
        }
        for obj in &self.marched_objects {
            let d = obj.sdf(p);//Slow vtable call
            if d < max_dis {
                max_dis = d;
                normal = obj.get_outward_normal(p);//Slow vtable call
                material = Some(*obj.material());//Slow vtable call
            }
        }
        return (max_dis,normal,material);
    }
}