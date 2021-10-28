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
}

pub struct HittableList{
    pub spheres: Vec<Sphere>,
    pub infinite_planes: Vec<InfinitePlane>,
    pub parallelograms: Vec<Parallelogram>,
    pub triangles: Vec<Triangle>,
    pub cubes: Vec<Cube>,
    pub objects: Vec<Box<dyn Traced + Send + Sync>>,
    pub marched_spheres: Vec<MarchedSphere>,
    pub marched_boxes: Vec<MarchedBox>,
    pub marched_torus: Vec<MarchedTorus>,
    pub marched_objects: Vec<Box<dyn Marched + Send + Sync>>,
}

macro_rules! add_obj {
    ( $x: ty, $name: ident ) => {
        impl std::ops::AddAssign<&$x> for HittableList{
            fn add_assign(&mut self, obj: &$x){
                self.$name.push(*obj)
            }
        }
    };
}

add_obj!(Sphere,spheres);
add_obj!(InfinitePlane,infinite_planes);
add_obj!(Parallelogram,parallelograms);
add_obj!(Triangle,triangles);
add_obj!(Cube,cubes);
add_obj!(MarchedSphere,marched_spheres);
add_obj!(MarchedBox,marched_boxes);
add_obj!(MarchedTorus,marched_torus);

macro_rules! new_and_clear {
    ($($name:ident),*) => {
        impl HittableList {
            pub fn new() -> Self{
                return HittableList{
                $(
                    $name: Vec::new(),
                )*
                };
            }
            #[allow(dead_code)]
            pub fn clear(&mut self) -> (){
                $(
                    self.$name.clear();
                )*
            }
        }
    };
}

new_and_clear!(spheres,infinite_planes,parallelograms,triangles,cubes,objects,marched_spheres,marched_boxes,marched_torus,marched_objects);

impl HittableList{
    #[allow(dead_code)]
    pub fn add_traced(&mut self,obj: Box<dyn Traced + Send + Sync>) -> () {
        self.objects.push(obj);
    }
    #[allow(dead_code)]
    pub fn add_marched(&mut self,obj: Box<dyn Marched + Send + Sync>) -> () {
        self.marched_objects.push(obj);
    }
    pub fn freeze(&self,cam: &Camera) -> FrozenHittableList{
        return FrozenHittableList::new(self,cam);
    }
}

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
        for obj in &self.hl.spheres{
            if let Some(hr) = obj.hit(r,t_min,closest_so_far) {
                closest_so_far = hr.t;
                rec = Some(hr);
            }
        }
        for obj in &self.hl.infinite_planes{
            if let Some(hr) = obj.hit(r,t_min,closest_so_far) {
                closest_so_far = hr.t;
                rec = Some(hr);
            }
        }
        for obj in &self.hl.parallelograms{
            if let Some(hr) = obj.hit(r,t_min,closest_so_far) {
                closest_so_far = hr.t;
                rec = Some(hr);
            }
        }
        for obj in &self.hl.triangles{
            if let Some(hr) = obj.hit(r,t_min,closest_so_far) {
                closest_so_far = hr.t;
                rec = Some(hr);
            }
        }
        for obj in &self.hl.cubes{
            if let Some(hr) = obj.hit(r,t_min,closest_so_far) {
                closest_so_far = hr.t;
                rec = Some(hr);
            }
        }
        for obj in &self.hl.objects{
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
            let (d,_,_,obj) = self.get_closest_distance_normal_material(&r.at(t));
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
            let (d,outward_normal,material,_) = self.get_closest_distance_normal_material(&p);

            //Should never happen the only raymarched object gets deleted mid transition between unstucking and raymarching
            if material.is_none() { return rec; }

            if  d < HIT_SIZE {//We hit something
                rec = Some(HitRecord{t: t,point: p,normal: outward_normal, material: material.unwrap()});
                break;
            }
            else { //Move forward
                t += d;//This only works if our direction in our Ray is unit length!!!
            }
        }
        return rec; 
    }

    pub fn get_closest_distance_normal_material(&self,p: &Point3) -> (f32,Vec3,Option<Material>,Option<Box<(dyn Marched + Send + Sync)>>){
        let mut max_dis    = INF;
        let mut normal = Vec3::ZERO;
        let mut material = None;
        let mut found_obj: Option<Box<(dyn Marched + Send + Sync)>> = None;
        for obj in &self.hl.marched_spheres {
            let d = obj.sdf(p).abs();//Not an actual vtable call, just a normal fast function call
            if d < max_dis {
                max_dis = d;
                normal = obj.get_outward_normal(p);//Not an actual vtable call, just a normal fast function call
                material = Some(obj.material);
                found_obj = Some(Box::new(*obj));
            }
        }
        for obj in &self.hl.marched_boxes {
            let d = obj.sdf(p).abs();//Not an actual vtable call, just a normal fast function call
            if d < max_dis {
                max_dis = d;
                normal = obj.get_outward_normal(p);//Not an actual vtable call, just a normal fast function call
                material = Some(obj.material);
                found_obj = Some(Box::new(*obj));
            }
        }
        for obj in &self.hl.marched_torus {
            let d = obj.sdf(p).abs();//Not an actual vtable call, just a normal fast function call
            if d < max_dis {
                max_dis = d;
                normal = obj.get_outward_normal(p);//Not an actual vtable call, just a normal fast function call
                material = Some(obj.material);
                found_obj = Some(Box::new(*obj));
            }
        }
        //@TODO: See how I can return dynamic Marched objects
        /*for obj in &self.marched_objects {
            let d = obj.sdf(p).abs();//Not an actual vtable call, just a normal fast function call
            if d < max_dis {
                max_dis = d;
                let o = obj.clone();
                normal = o.get_outward_normal(p);//Slow vtable call
                material = Some(*o.material());//Slow vtable call
                found_obj = Some(o);//I don't know if this is safe
            }
        }*/
        return (max_dis,normal,material,found_obj);
    }
}