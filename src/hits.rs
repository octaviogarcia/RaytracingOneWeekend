use crate::math::vec3::{Vec3,UnitVec3,Point3};
use crate::utils::{INF,get_id};
use crate::ray::Ray;
use crate::materials::Material;
use crate::traced::*;
use crate::marched::*;
use crate::camera::Camera;
use crate::bounding_box::Bounded;
use crate::camera_hash::CameraHash;

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

#[derive(Clone)]
pub struct HittableList{
    traced_objects: Vec<Arc<dyn Traced + Send + Sync>>,
    marched_objects: Vec<Arc<dyn Marched + Send + Sync>>,
    $($traced_ident: Vec<$traced>,)*
    $($marched_ident: Vec<$marched>,)*
}

pub struct FrozenHittableList{
    traced_objects: Vec<Arc<dyn Traced + Send + Sync>>,
    marched_objects: Vec<Arc<dyn Marched + Send + Sync>>,
    $($traced_ident: Vec<$traced>,)*
    $($marched_ident: Vec<$marched>,)*
    m_world_to_camera: crate::math::mat4x4::Mat4x4,
    #[allow(dead_code)]
    camera_hash: CameraHashes, 
}

pub struct CameraHashes{
    #[allow(dead_code)]
    traced_objects: CameraHash,
    #[allow(dead_code)]
    marched_objects: CameraHash,
    $(#[allow(dead_code)]
    $traced_ident: CameraHash,
    )*
    $(#[allow(dead_code)]
    $marched_ident: CameraHash,
    )*
}

impl CameraHashes{
    pub fn new() -> Self{
        Self {
            traced_objects: CameraHash::new(),
            marched_objects: CameraHash::new(),
            $($traced_ident: CameraHash::new(),)*
            $($marched_ident: CameraHash::new(),)*
        }
    }
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
    pub fn freeze(&mut self,cam: &Camera) -> FrozenHittableList{
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


const HIT_SIZE: f32 = 0.001;

impl FrozenHittableList{
    pub fn new(hl: &mut HittableList,cam: &Camera) -> Self{
        let viewmat = cam.viewmatrix();
        let viewmat_inv = viewmat.fast_homogenous_inverse();

        let mut ret = Self{
            traced_objects: hl.traced_objects.clone(),
            marched_objects: hl.marched_objects.clone(),
            $($traced_ident: hl.$traced_ident.clone(),)*
            $($marched_ident: hl.$marched_ident.clone(),)*
            m_world_to_camera: viewmat_inv,
            camera_hash: CameraHashes::new(),
        }; 
        $(
            for obj in &mut ret.$traced_ident{
                obj.build_bounding_box(&viewmat_inv);
            }
        )*
        for obj in &mut ret.traced_objects{//This modifies the base object rather than clone it, not sure how I feel about it
           Arc::get_mut(obj).unwrap().build_bounding_box(&viewmat_inv);
        }
        $(
            for obj in &mut ret.$marched_ident{
                obj.build_bounding_box(&viewmat_inv);
            }
        )*
        for obj in &mut ret.marched_objects {
            Arc::get_mut(obj).unwrap().build_bounding_box(&viewmat_inv);
        }
        return ret;
    }

    pub fn hit(&self,first_hit: bool,r: &Ray,t_min: f32,t_max: f32) -> Option<HitRecord> {
        //Ray tracing section
        let mut closest_so_far = t_max;
        let mut rec: Option<HitRecord>  = None;
        //If something isn't rendering properly, it might be because its not checking t_min,t_max in hit()
        let dir_from_camera = self.m_world_to_camera.dot_v3(&r.dir).to_z1();
        $(for obj in &self.$traced_ident{
            //Short circuit when not the first_hit, avoid calling hit_bounding_box
            let hit_bb = !first_hit || obj.hit_bounding_box(&dir_from_camera);
            if !hit_bb { continue; }
            if let Some(hr) = obj.hit(r,t_min,closest_so_far) {
                closest_so_far = hr.t;
                rec = Some(hr);
            }
        })*
        for obj in &self.traced_objects{
            let hit_bb = !first_hit || obj.hit_bounding_box(&dir_from_camera);
            if !hit_bb { continue; }
            if let Some(hr) = obj.hit(r,t_min,closest_so_far) {
                closest_so_far = hr.t;
                rec = Some(hr);
            }
        }

        //Ray marching section
        let mut t = self.unstuck(t_min,r);//If we started stuck in a object... unstuck ourselves
        if t.is_infinite() {//No marched objects in the scene. return raycasted result
            return rec;
        }
        let mut max_march_iter = 1024;

        //@SPEED: Replace this with a geohash at construction
        $(
            let mut $marched_ident: ([usize; 16],usize) = ([0;16],0);
            let mut idx = 0;
            for obj in &self.$marched_ident{
                let hit_bb = !first_hit || obj.hit_bounding_box(&dir_from_camera);
                if hit_bb {
                    $marched_ident.0[$marched_ident.1] = idx;
                    $marched_ident.1 += 1;
                }
                idx+=1;
            }
        )*
        let mut marched_objects: ([usize; 16],usize) = ([0;16],0);
        let mut idx = 0;
        for obj in &self.marched_objects{
            let hit_bb = !first_hit || obj.hit_bounding_box(&dir_from_camera);
            if hit_bb {
                marched_objects.0[marched_objects.1] = idx;
                marched_objects.1 += 1;
            }
            idx+=1;
        }

        while t < t_max && t < closest_so_far && max_march_iter > 0 {
            max_march_iter-=1;
            let point = r.at(t);
            let mut distance = INF;
            let mut normal   = Vec3::ZERO;
            let mut id       = 0;
            let mut material: Option<Material> = None;
            $(
                for idx_idx in 0..$marched_ident.1{
                    let idx = $marched_ident.0[idx_idx];
                    let d = self.$marched_ident[idx].sdf(&point).abs();//Not an actual vtable call, just a normal fast function call
                    if d < distance {
                        distance = d;
                        normal = self.$marched_ident[idx].get_outward_normal(&point);//Not an actual vtable call, just a normal fast function call
                        material = Some(self.$marched_ident[idx].material);
                        id = get_id(&self.$marched_ident[idx]);
                    }
                }
            )*
            for idx_idx in 0..marched_objects.1{
                let idx = marched_objects.0[idx_idx];
                let d = self.marched_objects[idx].sdf(&point).abs();//Not an actual vtable call, just a normal fast function call
                if d < distance {
                    distance = d;
                    normal = self.marched_objects[idx].get_outward_normal(&point);//Not an actual vtable call, just a normal fast function call
                    material = Some(*self.marched_objects[idx].material());
                    id = get_id(self.marched_objects[idx].as_ref());
                }
            }

            //Should never happen the only raymarched object gets deleted mid transition between unstucking and raymarching
            if material.is_none() { return rec; }

            if distance < HIT_SIZE {//We hit something
                rec = Some(HitRecord{t: t,point: point,normal: normal, material: material.unwrap(),obj_id: id});
                break;
            }
            else { //Move forward
                t += distance;//This only works if our direction in our Ray is unit length!!!
            }
        }
        return rec; 
    }

    pub fn unstuck(&self,t: f32,r: &Ray) -> f32{
        const MIN_STEP_SIZE: f32 = HIT_SIZE/2.;
        let mut new_t = t;
        let mut d = f32::INFINITY;
        let p = r.at(t);
        let mut marched: Option<&(dyn Marched + Send + Sync)> = None;

        $(for obj in &self.$marched_ident{
            let nd = obj.sdf(&p).abs();//Not an actual vtable call, just a normal fast function call
            if nd < d {
                d = nd;
                marched = Some(obj);
            }
        })*
        for obj in &self.marched_objects {
            let nd = obj.sdf(&p).abs();
            if nd < d {
                d = nd;
                marched = Some(obj.as_ref());
            }
        }

        let mut aux = d;
        if marched.is_none() { return f32::INFINITY; }//No marched objects on the Scene... return what we already
        while aux < HIT_SIZE {
            new_t += MIN_STEP_SIZE;
            aux = marched.unwrap().sdf(&r.at(new_t)).abs();
        }
        return new_t;
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
