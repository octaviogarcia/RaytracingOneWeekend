
use crate::math::vec3::{Vec3,UnitVec3,Point3};
use crate::math::mat4x4::Mat4x4;
use crate::ray::Ray;
use crate::materials::Material;
use crate::bounding_box::*;

pub trait Marched: Bounded {
    fn sdf(&self,p: &Point3) -> f32;
    fn get_outward_normal(&self,p: &Point3) -> Vec3;
    fn material(&self) -> &Material;
    fn center(&self) -> Point3;
}

pub fn get_outward_numeric_normal<M: Marched>(marched: &M,p: &Point3) -> UnitVec3{
    let eps = 0.0000001;
    let ex = Point3::new(eps, 0., 0.);
    let ey = Point3::new( 0.,eps, 0.);
    let ez = Point3::new( 0., 0.,eps);
    let x = marched.sdf(&(*p+ex)) - marched.sdf(&(*p-ex));
    let y = marched.sdf(&(*p+ey)) - marched.sdf(&(*p-ey));
    let z = marched.sdf(&(*p+ez)) - marched.sdf(&(*p-ez));
    let normal = Vec3::new(x,y,z).unit();

    //Flip the sign so always the SDF grows in the direction of the normal
    let test_ray = Ray::new(&marched.center(),&normal);
    let start     = test_ray.at(0.);
    let start_val = marched.sdf(&start);
    let end     = test_ray.at(1.);
    let end_val = marched.sdf(&end);
    let sign = [-1.,1.][(end_val > start_val) as usize];//If it grows, keep the sign. Else flip it
    return normal*sign;
}

#[derive(Copy, Clone)]
pub struct MarchedSphere {
    pub center: Point3,
    pub radius: f32,
    pub material: Material
}

impl Marched for MarchedSphere {
    fn sdf(&self,p: &Point3) -> f32 {
        return (*p - self.center).length() - self.radius;
    }
    fn get_outward_normal(&self,p: &Point3) -> Vec3 {
        let normal = (*p - self.center).unit();
        return normal;
    }
    fn material(&self) -> &Material{
        return &self.material;
    }
    fn center(&self) -> Point3{
        return self.center;
    }
}
impl Bounded for MarchedSphere {}

#[derive(Copy, Clone)]
pub struct MarchedBox {
    pub center: Point3,
    pub sizes: Vec3,
    pub material: Material
}

impl Marched for MarchedBox {
    fn sdf(&self,p: &Point3) -> f32 {
        let q = (&(*p-self.center)).abs() - self.sizes;
        return q.max(&Vec3::ZERO).length() + q.x().max(q.y().max(q.z())).min(0.);
    }
    fn get_outward_normal(&self,p: &Point3) -> Vec3 {
        return get_outward_numeric_normal(self,p);
    }
    fn material(&self) -> &Material{
        return &self.material;
    }
    fn center(&self) -> Point3{
        return self.center;
    }
}
impl Bounded for MarchedBox {}

#[derive(Copy, Clone)]
pub struct MarchedTorus {
    pub m_local_to_world_translate_rotate: Mat4x4,
    pub m_world_to_local_translate_rotate: Mat4x4,
    pub m_local_to_world_scale: Vec3,
    pub m_world_to_local_scale: Vec3,
    pub sizes: Vec3,//Vec2... actualy
    pub material: Material
}

impl MarchedTorus {
    pub fn new(m_local_to_world: &Mat4x4,local_sizes: &Vec3,mat: &Material) -> Self{
        let (t,r,s) = m_local_to_world.decompose_into_translate_rotate_scale();//TRS
        let translate_rotate = t.dot_mat(&r);//TR
        let scale = s.diag().xyz();
        let scale_inv = Vec3::new(1./scale.x(),1./scale.y(),1./scale.z());
        Self{
            m_local_to_world_translate_rotate: translate_rotate,//TR
            m_world_to_local_translate_rotate: translate_rotate.fast_homogenous_inverse(),//R^-1 T^-1
            m_local_to_world_scale: scale,
            m_world_to_local_scale: scale_inv,
            sizes: *local_sizes,
            material: *mat
        }
    }
    fn local_sdf(&self,p: &Point3) -> f32 {
        let p2 = *p;//- center, that is Point3::ZERO in local coords
        let q = Vec3::new(Vec3::new(p2.x(),p2.z(),0.).length()-self.sizes.x(),p2.y(),0.);
        return q.length() - self.sizes.y();
    }
}

impl Marched for MarchedTorus {
    fn sdf(&self,p: &Point3) -> f32 {
        //http://jamie-wong.com/2016/07/15/ray-marching-signed-distance-functions/#non-uniform-scaling-and-beyond
        let local_point = self.m_world_to_local_translate_rotate.dot_p3(&(*p*self.m_world_to_local_scale));
        let min_scale = self.m_local_to_world_scale.min_val();//@SPEED: if non uniform scaling, the marching is suboptimal
        return self.local_sdf(&local_point)*min_scale;
    }
    fn get_outward_normal(&self,p: &Point3) -> Vec3 {
        return get_outward_numeric_normal(self,p);
    }
    fn material(&self) -> &Material{
        return &self.material;
    }
    fn center(&self) -> Point3{
        return self.m_local_to_world_translate_rotate.dot_p3(&Point3::ZERO);
    }
}
impl Bounded for MarchedTorus {}