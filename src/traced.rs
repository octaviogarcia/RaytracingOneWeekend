

use crate::vec3::{Vec3,UnitVec3,Point3};
use crate::mat3x3::{Mat3x3};
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

#[derive(Copy, Clone)]
pub struct Barycentric<const BT: usize>{
    pub origin: Point3,
    pub u: UnitVec3,
    pub u_length: f32,
    pub v: UnitVec3,
    pub v_length: f32,
    pub uxv: UnitVec3,
    pub base_inv: Mat3x3,//Transforms a point (X,Y,Z) to the plane 2d coordinates
    pub v_in_base: Vec3,//Vec2, Z is 0 @SPEED
    pub material: Material,
}

//Current rust version doesn't suport direct const enum templates... so I use an usize...
impl <const BT: usize> Barycentric<BT> {
    pub fn new(origin: &Point3,u: &UnitVec3,v: &UnitVec3,u_length: f32,v_length: f32,material: &Material) -> Self{
        let u_unit = u.unit();
        let v_unit = v.unit();
        let uxv = u_unit.cross(v_unit).unit();
        let uxvxu = uxv.cross(u_unit).unit();
        let base_inv = Mat3x3::new3v_vert(&u_unit,&uxv,&uxvxu).inverse();
        let v_in_base = base_inv.dot(&v_unit);
        return Self{origin: *origin,material: *material,
            u: u_unit,u_length: u_length,v: v_unit,v_length: v_length,uxv: uxv,
            base_inv: base_inv, v_in_base: v_in_base
        };
    }
    pub fn new3points(origin: &Point3,upoint: &Point3,vpoint: &Point3,material: &Material) -> Self{
        let u_rel = *upoint-*origin;
        let u_length = u_rel.length();
        let v_rel = *vpoint-*origin;
        let v_length = v_rel.length();
        return Self::new(origin,&u_rel,&v_rel,u_length,v_length,material);
    }
    #[inline]
    fn calc_barycentric(&self,coords2d: &Point3) -> (f32,f32,f32){
        let rx = coords2d.x();
        let ry = coords2d.z();
        let ux = 1.;
        let uy = 0.;
        let vx = self.v_in_base.x();
        let vy = self.v_in_base.z();
        let det = ux*vy - vx*uy;
        let lambda1 = ( (rx*vy - vx*ry)/det)/self.u_length;
        let lambda2 = (-(rx*uy - ux*ry)/det)/self.v_length;
        return (lambda1,lambda2,1.-lambda1-lambda2);
    }
    #[inline]
    fn check_lambdas_triangle(&self,l1: f32,l2: f32,l3: f32) -> bool{
        return l1 > 0. && l2 > 0. && l3 > 0. && l1 < 1. && l2 < 1. && l3 < 1.;
    }
    #[inline]
    fn check_lambdas_parallelogram(&self,l1: f32,l2: f32,l3: f32) -> bool{
        return l1 > 0. && l2 > 0. && l1 < 1. && l2 < 1.;
    }
    fn hit_aux(&self,r: &Ray,t_min: f32,t_max: f32) -> Option<HitRecord> {
        let (root,normal_dot_dir) = ray_plane_intersect(r,&self.uxv,&self.origin);
        if root == INF || root < t_min || root > t_max {
            return None;
        }
        let point = r.at(root);
        let point_from_origin = point - self.origin;
        let coords2d = self.base_inv.dot(&point_from_origin);//@SPEED: Y coord is useless, its 0 since we intersected
        let (lambda1, lambda2, lambda3) = self.calc_barycentric(&coords2d);
        let correct: bool;
        if BT == 0 {//Hope the compiler optimizes these static IFs!
            correct = self.check_lambdas_parallelogram(lambda1,lambda2,lambda3);
        }
        else if BT == 1 {
            correct = self.check_lambdas_triangle(lambda1,lambda2,lambda3);
        }
        else{
            panic!("Barycentric case not fulfilled");
        }
        if !correct { 
            return None;
        }
        let outward_normal = normal_against_direction(&self.uxv,normal_dot_dir);
        return Some(HitRecord{t: root,point: point,normal: outward_normal,material: self.material});
    }
}

impl <const BT: usize> Traced for Barycentric<BT> {
    fn hit(&self,r: &Ray,t_min: f32,t_max: f32) -> Option<HitRecord> {
        return self.hit_aux(r,t_min,t_max);
    }
}

pub type Parallelogram = Barycentric<0>;
pub type Triangle = Barycentric<1>;