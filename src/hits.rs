use crate::math::vec3::{Vec3,UnitVec3,Point3};
use crate::utils::INF;
use crate::ray::Ray;
use crate::materials::Material;
use crate::traced::*;
use crate::marched::*;
use crate::camera::Camera;

pub struct HitRecord {
    pub point: Point3,
    pub normal: UnitVec3,//Always outward from the surface
    pub material: Material,
    pub t: f32,
    pub obj_id: u64,
}

use std::sync::Arc;

macro_rules! hittable_list {
($($traced_ident:ident ; $traced:ty),* | $($marched_ident:ident ; $marched:ty),*) => {

pub struct HittableList{
    pub traced_objects: Vec<Arc<dyn Traced + Send + Sync>>,
    pub marched_objects: Vec<Arc<dyn Marched + Send + Sync>>,
    $($traced_ident: Vec<$traced>,)*
    $($marched_ident: Vec<$marched>,)*
}
impl HittableList {
    pub fn new() -> Self{
        return HittableList{
            traced_objects: Vec::new(),
            marched_objects: Vec::new(),
            $($traced_ident: Vec::new(),)*
            $($marched_ident: Vec::new(),)*
        };
    }
    #[allow(dead_code)]
    pub fn clear(&mut self) -> (){
        self.traced_objects.clear();
        self.marched_objects.clear();
        $(self.$traced_ident.clear();)*
        $(self.$marched_ident.clear();)*
    }
    pub fn freeze(&self,cam: &Camera) -> FrozenHittableList{
        return FrozenHittableList::new(self,cam);
    }
}
impl std::ops::AddAssign<Arc<dyn Traced + Send + Sync>> for HittableList {
    fn add_assign(&mut self, obj: Arc<dyn Traced + Send + Sync>){
        self.traced_objects.push(obj);
    }
}
impl std::ops::AddAssign<Arc<dyn Marched + Send + Sync>> for HittableList {
    fn add_assign(&mut self, obj: Arc<dyn Marched + Send + Sync>){
        self.marched_objects.push(obj);
    }
}
$(impl std::ops::AddAssign<&$traced> for HittableList{
    fn add_assign(&mut self, obj: &$traced){
        self.$traced_ident.push(*obj);
    }
})*
$(impl std::ops::AddAssign<&$marched> for HittableList{
    fn add_assign(&mut self, obj: &$marched){
        self.$marched_ident.push(*obj);
    }
})*

pub struct FrozenHittableList<'a>{
    hl: &'a HittableList,//@TODO: implement some sort of KD tree or octotree
}

impl <'a> FrozenHittableList<'a>{
    pub fn new(hl: &'a HittableList,_cam: &Camera) -> Self{ Self{hl: hl} }
    #[allow(dead_code)]
    pub fn unfreeze(&self) -> &HittableList { self.hl }
    pub fn hit(&self,r: &Ray,t_min: f32,t_max: f32) -> Option<HitRecord> {
        //Ray tracing section
        let mut closest_so_far = t_max;
        let mut rec: Option<HitRecord>  = None;
        //If something isn't rendering properly, it might be because its not checking t_min,t_max in hit()
        $(for obj in &self.hl.$traced_ident{
            if let Some(hr) = obj.hit(r,t_min,closest_so_far) {
                closest_so_far = hr.t;
                rec = Some(hr);
            }
        })*
        for obj in &self.hl.traced_objects{
            if let Some(hr) = obj.hit(r,t_min,closest_so_far) {
                closest_so_far = hr.t;
                rec = Some(hr);
            }
        }

        //Ray marching section
        const HIT_SIZE: f32 = 0.001;
        let mut t = t_min;
        {//If we started stuck in a wall... unstuck ourselves
            const MIN_STEP_SIZE: f32 = HIT_SIZE/2.;
            let (d,_,_,obj,_) = self.get_closest_distance_normal_material(&r.at(t));
            let mut aux = d;
            if obj.is_none() { return rec; }
            let o = obj.unwrap();//Maybe its faster to match the Option than unwrap()... I'm not sure
            while aux < HIT_SIZE {
                t += MIN_STEP_SIZE;
                aux = o.sdf(&r.at(t)).abs();
            }
        }

        let mut max_march_iter = 1024;
        while t < t_max && t < closest_so_far && max_march_iter > 0 {
            max_march_iter-=1;
            let p = r.at(t);
            let (d,outward_normal,material,_,obj_id) = self.get_closest_distance_normal_material(&p);

            //Should never happen the only raymarched object gets deleted mid transition between unstucking and raymarching
            if material.is_none() { return rec; }

            if  d < HIT_SIZE {//We hit something
                //let obj_id: u64 = std::ptr::addr_of!(obj) as u64;
                rec = Some(HitRecord{t: t,point: p,normal: outward_normal, material: material.unwrap(),obj_id});
                break;
            }
            else { //Move forward
                t += d;//This only works if our direction in our Ray is unit length!!!
            }
        }
        return rec; 
    }

    pub fn get_closest_distance_normal_material(&self,p: &Point3) -> (f32,Vec3,Option<Material>,Option<Arc<(dyn Marched + Send + Sync)>>,u64){
        let mut max_dis    = INF;
        let mut normal = Vec3::ZERO;
        let mut material = None;
        let mut found_obj: Option<Arc<(dyn Marched + Send + Sync)>> = None;
        let mut id: u64 = 0;
        $(for obj in &self.hl.$marched_ident{
            let d = obj.sdf(p).abs();//Not an actual vtable call, just a normal fast function call
            if d < max_dis {
                max_dis = d;
                normal = obj.get_outward_normal(p);//Not an actual vtable call, just a normal fast function call
                material = Some(obj.material);
                found_obj = Some(Arc::new(*obj));//@SPEED This is a clone...
                id = obj.get_id();
            }
        })*
        for obj in &self.hl.marched_objects {
            let d = obj.sdf(p).abs();
            if d < max_dis {
                max_dis = d;
                normal = obj.get_outward_normal(p);//Slow vtable call
                material = Some(*obj.material());//Slow vtable call
                found_obj = Some(obj.clone());
                id = obj.get_id();
            }
        }
        return (max_dis,normal,material,found_obj,id);
    }
}
 
};}

hittable_list!(spheres;Sphere, cubes;Cube, triangles;Triangle,infinite_planes;InfinitePlane,parallelograms;Parallelogram
              |marched_spheres;MarchedSphere,marched_boxes;MarchedBox,marched_torus;MarchedTorus);

//Only tested with abs(obj.sdf(r.at(t))) < HIT_SIZE
#[allow(dead_code)]
fn root_find(obj: Option<Box<(dyn Marched + Send + Sync)>>,t: f32,r: &Ray,hit_size: f32) -> f32 {
    let o = obj.unwrap();

    let mut first_side_t   = t;
    let mut first_side_val = o.sdf(&r.at(t));
    //If positive ADD to the ray, (we are going inside the surface). Else, SUBSTRACT, we are going outside.
    let sign = [-1.,1.][(first_side_val > 0.) as usize];
    let mut other_side_t   = first_side_t + (sign)*hit_size;
    let mut other_side_val = o.sdf(&r.at(other_side_t));
    let mut max_iters = 10;
    while (first_side_val > 0.) == (other_side_val > 0.)  && max_iters > 0{//Find a point on the other side
        other_side_t += (sign)*hit_size;
        other_side_val = o.sdf(&r.at(other_side_t));
        max_iters-=1;
    }
    if max_iters == 0 {return INF;}
    //Now bisect
    let mut middle_t   = (first_side_t + other_side_t)/2.;
    let mut middle_val = o.sdf(&r.at(middle_t));
    for _i in 0..5 {
        if (middle_val > 0.) == (first_side_val > 0.){
            first_side_t   = middle_t;
            first_side_val = middle_val;
        }
        else if (middle_val > 0.) == (other_side_val > 0.){
            other_side_t   = middle_t;
            other_side_val = middle_val;
        }
        middle_t   = (first_side_t + other_side_t)/2.;
        middle_val = o.sdf(&r.at(middle_t));
    }
    //Asume its good enough to do Newtons
    let mut ti = middle_t;
    let mut fi = middle_val;
    for _i in 0..20 {
        let eps = hit_size/10.;
        let feps = o.sdf(&r.at(ti+eps));
        let dfi = (feps - fi)/eps;
        ti = ti - fi/dfi;
        fi = o.sdf(&r.at(ti)); 
    }
    return ti;
}
