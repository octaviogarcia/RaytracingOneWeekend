use crate::math::vec3::{Vec3,UnitVec3,Point3};
use crate::utils::{INF,get_id,VecIndexes};
use crate::ray::Ray;
use crate::materials::Material;
use crate::traced::*;
use crate::marched::*;
use crate::camera::Camera;
use crate::bounding_box::Bounded;
use crate::camera_hash::*;

pub struct HitRecord {
    pub point: Point3,
    pub normal: UnitVec3,//Always outward from the surface
    pub material: Material,
    pub t: f32,
    pub obj_id: u64,
}

use std::sync::Arc;


macro_rules! set_hash_index {
    ($camera_hash:expr,$bb:expr,$idx:expr,$ident:ident) => {
        let (mini,maxi,minj,maxj) = $camera_hash.get_indexes(&$bb);
        for i in mini..=maxi{
            for j in minj..=maxj{
                $camera_hash.cells[i][j].$ident.add($idx);
            }
        }
    }
}

macro_rules! hittable_list {
($($traced_ident:ident ; $traced:ty),* | $($marched_ident:ident ; $marched:ty),*) => {

#[derive(Clone)]
pub struct HittableList{
    traced_objects: Vec<Arc<dyn Traced + Send + Sync>>,
    marched_objects: Vec<Arc<dyn Marched + Send + Sync>>,
    $($traced_ident: Vec<$traced>,)*
    $($marched_ident: Vec<$marched>,)*
}

const MAX_INDEXES_HASH: usize = 16;
#[derive(Copy,Clone,Debug)]
struct CameraHashCell {
    pub traced_objects: VecIndexes<MAX_INDEXES_HASH>,
    pub marched_objects: VecIndexes<MAX_INDEXES_HASH>,
    $(pub $traced_ident: VecIndexes<MAX_INDEXES_HASH>,)*
    $(pub $marched_ident: VecIndexes<MAX_INDEXES_HASH>,)*
}

impl CellEmptyInitializable for CameraHashCell {
    fn empty(&mut self) -> () {
        self.traced_objects.empty();
        self.marched_objects.empty();
        $(self.$traced_ident.empty();)*
        $(self.$marched_ident.empty();)*
    }
}

pub struct FrozenHittableList{
    traced_objects: Vec<Arc<dyn Traced + Send + Sync>>,
    marched_objects: Vec<Arc<dyn Marched + Send + Sync>>,
    $($traced_ident: Vec<$traced>,)*
    $($marched_ident: Vec<$marched>,)*
    camera_hash: Box<CameraHash<CameraHashCell>>, 
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
        let world_to_camera = cam.world_to_camera();

        let mut ret = Self{
            traced_objects: hl.traced_objects.clone(),
            marched_objects: hl.marched_objects.clone(),
            $($traced_ident: hl.$traced_ident.clone(),)*
            $($marched_ident: hl.$marched_ident.clone(),)*
            camera_hash: unsafe {//Rust is just awful, all this mess just to init on the heap
                let layout = std::alloc::Layout::new::<CameraHash<CameraHashCell>>();
                let ptr = std::alloc::alloc(layout) as *mut CameraHash<CameraHashCell>;
                (*ptr).empty(); 
                Box::from_raw(ptr)
            },
        };

        //ret.camera_hash.bb 
        let bb = cam.viewport_world_bounding_box();
        let bb = bb.dot(&world_to_camera);
        let bb = bb.project(cam.focus_dist);//[-viewport*0.5;viewport*0.5]
        let bb = bb.scale(2./cam.viewport_width,2./cam.viewport_height);//[-1;1]
        let bb = bb.translate(1.,1.);//[0;2]
        let bb = bb.scale(0.5,0.5);//--> should be {0.,0.,1.,1.}
        println!("{:?}",bb);
        $({
            let mut idx = 0;
            for obj in &mut ret.$traced_ident{
                let bb = obj.build_world_bounding_box();
                let bb = bb.dot(&world_to_camera);
                let bb = bb.project(cam.focus_dist);
                let bb = bb.scale(1./cam.viewport_width,1./cam.viewport_height);//[-1;1]
                let bb = bb.translate(1.,1.);//[0;2]
                let bb = bb.scale(0.5,0.5);
                println!("{:?}",bb);
                set_hash_index!(ret.camera_hash,bb,idx,$traced_ident);
                idx += 1;
            }
        })*

        {
            let mut idx = 0;
            for obj in &mut ret.traced_objects{//This modifies the base object rather than clone it, not sure how I feel about it
                let bb = Arc::get_mut(obj).unwrap().build_world_bounding_box();
                let bb = bb.dot(&world_to_camera);
                let bb = bb.project(cam.focus_dist);
                let bb = bb.scale(1./(2.*cam.viewport_width),1./(2.*cam.viewport_height));//[-1;1]
                let bb = bb.translate(1.,1.);//[0;2]
                let bb = bb.scale(0.5,0.5);
                set_hash_index!(ret.camera_hash,bb,idx,traced_objects);
                idx += 1;
            }
        }

        $({
            let mut idx = 0;
            for obj in &mut ret.$marched_ident{
                let bb = obj.build_world_bounding_box();
                let bb = bb.dot(&world_to_camera);
                let bb = bb.project(cam.focus_dist);
                let bb = bb.scale(1./(2.*cam.viewport_width),1./(2.*cam.viewport_height));//[-1;1]
                let bb = bb.translate(1.,1.);//[0;2]
                let bb = bb.scale(0.5,0.5);
                set_hash_index!(ret.camera_hash,bb,idx,traced_objects);
                idx += 1;
            }
        })*

        {
            let mut idx = 0;
            for obj in &mut ret.marched_objects {
                let bb = Arc::get_mut(obj).unwrap().build_world_bounding_box();
                let bb = bb.dot(&world_to_camera);
                let bb = bb.project(cam.focus_dist);
                let bb = bb.scale(1./(2.*cam.viewport_width),1./(2.*cam.viewport_height));//[-1;1]
                let bb = bb.translate(1.,1.);//[0;2]
                let bb = bb.scale(0.5,0.5);
                set_hash_index!(ret.camera_hash,bb,idx,traced_objects);
                idx += 1;
            }
        }

        return ret;
    }

    pub fn first_hit(&self,r: &Ray,t_min: f32,t_max: f32,u: f32,v: f32) -> Option<HitRecord> {
        //Ray tracing section
        let mut closest_so_far = t_max;
        let mut rec: Option<HitRecord>  = None;
        //println!("{},{} -> {},{}",u,v,r.dir.to_z1().x(),r.dir.to_z1().y());
        let cell = self.camera_hash.at(u,v);
        $({
            let arr = &cell.$traced_ident;
            for idx_idx in 0..arr.count{
                let idx = arr.arr[idx_idx];
                if let Some(hr) = self.$traced_ident[idx].hit(r,t_min,closest_so_far) {
                    closest_so_far = hr.t;
                    rec = Some(hr);
                }
            }
        })*
        {
            let arr = &cell.traced_objects;
            for idx_idx in 0..arr.count{
                let idx = arr.arr[idx_idx];
                if let Some(hr) = self.traced_objects[idx].hit(r,t_min,closest_so_far) {
                    closest_so_far = hr.t;
                    rec = Some(hr);
                }
            }
        }

        let mut t = self.unstuck(t_min,r);//If we started stuck in a object... unstuck ourselves
        if t.is_infinite() {//No marched objects in the scene. return raycasted result
            return rec;
        }
        let mut max_march_iter = 1024;

        while t < t_max && t < closest_so_far && max_march_iter > 0 {
            max_march_iter-=1;
            let point = r.at(t);
            let mut distance = INF;
            let mut normal   = Vec3::ZERO;
            let mut id       = 0;
            let mut material: Option<Material> = None;
            $({
                let arr = cell.$marched_ident;
                for idx_idx in 0..arr.count{
                    let idx = arr.arr[idx_idx];
                    let obj = &self.$marched_ident[idx];
                    let d = obj.sdf(&point).abs();//Not an actual vtable call, just a normal fast function call
                    if d < distance {
                        distance = d;
                        normal = obj.get_outward_normal(&point);//Not an actual vtable call, just a normal fast function call
                        material = Some(obj.material);
                        id = get_id(obj);
                    }
                }
            })*
            {
                let arr = cell.marched_objects;
                for idx_idx in 0..arr.count{
                    let idx = arr.arr[idx_idx];
                    let obj = &self.marched_objects[idx];
                    let d = obj.sdf(&point).abs();//Not an actual vtable call, just a normal fast function call
                    if d < distance {
                        distance = d;
                        normal = obj.get_outward_normal(&point);//Not an actual vtable call, just a normal fast function call
                        material = Some(*obj.material());
                        id = get_id(obj.as_ref());
                    }
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

    pub fn hit(&self,r: &Ray,t_min: f32,t_max: f32) -> Option<HitRecord> {
        //Ray tracing section
        let mut closest_so_far = t_max;
        let mut rec: Option<HitRecord>  = None;
        $(for obj in &self.$traced_ident{
            if let Some(hr) = obj.hit(r,t_min,closest_so_far) {
                closest_so_far = hr.t;
                rec = Some(hr);
            }
        })*
        for obj in &self.traced_objects{
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

        while t < t_max && t < closest_so_far && max_march_iter > 0 {
            max_march_iter-=1;
            let point = r.at(t);
            let mut distance = INF;
            let mut normal   = Vec3::ZERO;
            let mut id       = 0;
            let mut material: Option<Material> = None;
            $(
                for obj in &self.$marched_ident{
                    let d = obj.sdf(&point).abs();//Not an actual vtable call, just a normal fast function call
                    if d < distance {
                        distance = d;
                        normal = obj.get_outward_normal(&point);//Not an actual vtable call, just a normal fast function call
                        material = Some(obj.material);
                        id = get_id(obj);
                    }
                }
            )*
            for obj in &self.marched_objects{
                let d = obj.sdf(&point).abs();
                if d < distance {
                    distance = d;
                    normal = obj.get_outward_normal(&point);
                    material = Some(*obj.material());
                    id = get_id(obj.as_ref());
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
