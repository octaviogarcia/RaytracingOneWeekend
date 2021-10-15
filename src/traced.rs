

use crate::vec3;
use vec3::{Vec3,UnitVec3,Point3,Mat3x3};
use crate::utils::{INF};
use crate::ray::Ray;
use crate::hits::HitRecord;
use crate::materials::Material;

#[derive(Copy, Clone)]
pub struct Sphere {
    pub center: Point3,
    pub radius: f32,
    pub material: Material
}

pub trait Traced {
    fn hit(&self,r: &Ray,t_min: f32,t_max: f32) -> Option<HitRecord>;
}

impl Traced for Sphere {
    fn hit(&self,r: &Ray,t_min: f32,t_max: f32) -> Option<HitRecord> {
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
        return Some(HitRecord{t: root,point: point,normal: outward_normal,material: self.material});
    }
}

#[derive(Copy, Clone)]
pub struct InfinitePlane {
    pub center: Point3,
    pub normal: UnitVec3,
    pub material: Material,
}

impl InfinitePlane {
    pub fn new(center: &Point3,normal: &UnitVec3,material: &Material) -> Self{
        return InfinitePlane{center: *center,normal: normal.unit(),material: *material};
    }
}
/*
//Bunch of ifs to get expected values for the axis directions... probably a more consistent way...
if      normal.y() == 0. && normal.z() == 0.{//(nx 0 0) generates (0 nx 0),(0 0 nx)
    u = Vec3::new(0.,normal.x(),0.).unit();
    v = Vec3::new(0.,       0.,normal.x()).unit();
}
else if normal.x() == 0. && normal.z() == 0. {//(0 ny 0) generates (ny 0 0),(0 0 ny)
    u = Vec3::new(normal.y(),0.,0.).unit();
    v = Vec3::new(       0., 0.,normal.y()).unit();
}
else if normal.x() == 0. && normal.y() == 0. {//(0 0 nz) generates (nz 0 0),(0 nz 0)
    u = Vec3::new(normal.z(),       0.,0.).unit();
    v = Vec3::new(       0.,normal.z(),0.).unit();
}
else{
    u = Vec3::new(-normal.y(),normal.x(),0.).unit();
    v = normal.cross(u).unit();
}*/

#[inline]
fn ray_plane_intersect(r: &Ray,plane_normal: &UnitVec3,plane_center: &Point3) -> (f32,f32){
    let div = plane_normal.dot(r.dir);
    if div.abs() < 0.000001{//Plane is paralel to the ray
        return (INF,0.);
    }
    let num = -plane_normal.dot(r.orig - *plane_center);
    return (num/div,div);
}
#[inline]
fn normal_against_direction(normal: &Vec3,normal_dot_dir: f32) -> Vec3{
    return *normal*([1.,-1.][(normal_dot_dir > 0.) as usize]);//if positive is same direction... reverse it
}

impl Traced for InfinitePlane {
    fn hit(&self,r: &Ray,t_min: f32,t_max: f32) -> Option<HitRecord> {
        let (root,normal_dot_dir) = ray_plane_intersect(r,&self.normal,&self.center);
        if root == INF || root < t_min || root > t_max {
            return None;
        }
        let outward_normal = normal_against_direction(&self.normal,normal_dot_dir);
        //Maybe its faster to send some sort of reference/pointer to material? Probably not, since its so small
        return Some(HitRecord{t: root,point: r.at(root),normal: outward_normal,material: self.material});
    }
}

pub type Triangle = Parallelogram;
#[derive(Copy, Clone)]
pub struct Parallelogram {//Easier in parametric [f(t,s) = u*t+v*s+origin] form
    pub origin: Point3,
    pub u: UnitVec3,
    pub u_length: f32,
    pub v: UnitVec3,
    pub v_length: f32,
    pub uxv: UnitVec3,
    pub ubase_inv: Mat3x3,//Transforms a point (X,Y,Z) to the plane local coordinates (UCOORD,_,_)
    pub vbase_inv: Mat3x3,//Transforms a point (X,Y,Z) to the plane local coordinates (_,VCOORD,_)
    pub material: Material,
}

impl Parallelogram {
    pub fn new(origin: &Point3,u: &UnitVec3,v: &UnitVec3,u_length: f32,v_length: f32,material: &Material) -> Self{
        let u_unit = u.unit();
        let v_unit = v.unit();
        let uxv = u_unit.cross(v_unit).unit();
        let ubase_inv = Mat3x3::new3v(&u_unit,&uxv,&uxv.cross(u_unit).unit()).inverse();
        let vbase_inv = Mat3x3::new3v(&v_unit,&uxv,&uxv.cross(v_unit).unit()).inverse();
        return Self{origin: *origin,material: *material,
            u: u_unit,u_length: u_length,v: v_unit,v_length: v_length,uxv: uxv,
            ubase_inv: ubase_inv,vbase_inv: vbase_inv};
    }
    pub fn new3points(origin: &Point3,upoint: &Point3,vpoint: &Point3,material: &Material) -> Self{
        let u_rel = *upoint-*origin;
        let u_length = u_rel.length();
        let u_unit = u_rel/u_length;

        let v_rel = *vpoint-*origin;
        let v_length = v_rel.length();
        let v_unit = v_rel/v_length;
        return Self::new(origin,&u_unit,&v_unit,u_length,v_length,material);
    }
}

impl Traced for Parallelogram {
    fn hit(&self,r: &Ray,t_min: f32,t_max: f32) -> Option<HitRecord> {
        let (root,normal_dot_dir) = ray_plane_intersect(r,&self.uxv,&self.origin);
        if root == INF || root < t_min || root > t_max {
            return None;
        }
        let point = r.at(root);
        let point_from_origin = point - self.origin;
        let ucoord = self.ubase_inv.dot(&point_from_origin).e[0];
        let vcoord = self.vbase_inv.dot(&point_from_origin).e[0];
        if ucoord < 0. || vcoord < 0. || ucoord > self.u_length || vcoord > self.v_length {
            return None;
        } 
        let outward_normal = normal_against_direction(&self.uxv,normal_dot_dir);
        return Some(HitRecord{t: root,point: point,normal: outward_normal,material: self.material});
    }
}

/*
#[derive(Copy, Clone)]
pub struct Cube {
    pub center: Point3,
    pub normals: [UnitVec3;6],
    pub length: f32,
    pub material: Material
}

impl Cube{
    pub fn new(center: &Point3,length: f32,material: &Material) -> Self {
        let n = {}
        return Self{center: *center,length: length,material: *material};
    }
}

impl Traced for Cube {
    fn hit(&self,r: &Ray,t_min: f32,t_max: f32) -> Option<HitRecord> {
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
        return Some(HitRecord{t: root,point: point,normal: outward_normal,material: self.material});
    }
}*/