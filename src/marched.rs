
use crate::math::vec3::{Vec3,UnitVec3,Point3};
use crate::math::mat4x4::Mat4x4;
use crate::math::vec4::Vec4;
use crate::ray::Ray;
use crate::materials::Material;
use crate::bounding_box::*;
pub trait Marched: Bounded {
    fn material(&self) -> &Material;
    fn to_local(&self,p: &Vec4) -> Vec4;
    fn to_world(&self,p: &Vec4) -> Vec4;
    fn to_world_f(&self,f: f32) -> f32;
    fn local_sdf(&self,p: &Point3) -> f32;
    fn sdf(&self,p: &Point3) -> f32{
        let p = self.to_local(&Vec4::new_p3(p)).xyz();
        let local_sdf = self.local_sdf(&p);
        return self.to_world_f(local_sdf);
    }
    fn get_outward_normal(&self,p: &Point3) -> UnitVec3{
        let p = self.to_local(&Vec4::new_p3(p)).xyz();
        let n = Vec4::new_v3(&self.get_outward_local_normal(&p));
        let world_n = self.to_world(&n);
        return world_n.xyz().unit();
    }
    fn get_outward_local_normal(&self,p: &Point3) -> UnitVec3 {
        let eps = 0.0000001;
        let ex = Point3::new(eps, 0., 0.);
        let ey = Point3::new( 0.,eps, 0.);
        let ez = Point3::new( 0., 0.,eps);
        let x = self.local_sdf(&(*p+ex)) - self.local_sdf(&(*p-ex));
        let y = self.local_sdf(&(*p+ey)) - self.local_sdf(&(*p-ey));
        let z = self.local_sdf(&(*p+ez)) - self.local_sdf(&(*p-ez));
        let normal = Vec3::new(x,y,z).unit();
        //Flip the sign so always the SDF grows in the direction of the normal
        // This should work in local and word coords as long its convex
        // and there is no inversion or shearing or something like that
        let test_ray = Ray::new(&Point3::ZERO,&normal);//Ray::new(&self.center(),&normal);
        let start     = test_ray.at(0.);
        let start_val = self.local_sdf(&start);
        let end     = test_ray.at(1.);
        let end_val = self.local_sdf(&end);
        let sign = [-1.,1.][(end_val > start_val) as usize];//If it grows, keep the sign. Else flip it
        return normal*sign;
    }
}
//For now, just always draw the marched


#[derive(Copy, Clone)]
pub struct MarchedSphere {
    pub center: Point3,
    pub radius: f32,
    pub material: Material
}
impl Bounded for MarchedSphere {}
impl Marched for MarchedSphere {
    fn local_sdf(&self,p: &Point3) -> f32 {
        return p.length() - self.radius;
    }
    fn get_outward_normal(&self,p: &Point3) -> Vec3 {
        let normal = (*p - self.center).unit();
        return normal;
    }
    fn material(&self) -> &Material{
        return &self.material;
    }
    fn to_local(&self,p: &Vec4) -> Vec4 {
        return *p - p.w()*Vec4::new_v3(&self.center);
    }
    fn to_world(&self,p: &Vec4) -> Vec4{
        return *p + p.w()*Vec4::new_v3(&self.center);
    }
    fn to_world_f(&self,f: f32) -> f32{
        return f;
    }
}

#[derive(Copy, Clone)]
pub struct MarchedBox {
    pub center: Point3,
    pub sizes: Vec3,
    pub material: Material
}
impl Bounded for MarchedBox {}
impl Marched for MarchedBox {
    fn local_sdf(&self,p: &Point3) -> f32 {
        let q = p.abs() - self.sizes;
        return q.max(&Vec3::ZERO).length() + q.x().max(q.y().max(q.z())).min(0.);
    }
    fn material(&self) -> &Material{
        return &self.material;
    }
    fn to_local(&self,p: &Vec4) -> Vec4 {
        return *p - p.w()*Vec4::new_v3(&self.center);
    }
    fn to_world(&self,p: &Vec4) -> Vec4{
        return *p + p.w()*Vec4::new_v3(&self.center);
    }
    fn to_world_f(&self,f: f32) -> f32{
        return f;
    }
}

#[derive(Copy, Clone)]
pub struct MarchedTorus {
    pub m_local_to_world_translate_rotate: Mat4x4,
    pub m_world_to_local_translate_rotate: Mat4x4,
    pub m_local_to_world_scale: Vec4,
    pub m_world_to_local_scale: Vec4,
    pub sizes: Vec3,//Vec2... actualy
    pub material: Material,
    pub bounding_box: BoundingBox,
}

impl MarchedTorus {
    pub fn new(m_local_to_world: &Mat4x4,local_sizes: &Vec3,mat: &Material) -> Self{
        let (t,r,s) = m_local_to_world.decompose_into_translate_rotate_scale();//TRS
        let translate_rotate = t.dot_mat(&r);//TR
        let scale = s.diag();//w value should be 1.
        let scale_inv = Vec4::new(1./scale.x(),1./scale.y(),1./scale.z(),1./scale.w());//w value should be 1.
        Self{
            m_local_to_world_translate_rotate: translate_rotate,//TR
            m_world_to_local_translate_rotate: translate_rotate.fast_homogenous_inverse(),//R^-1 T^-1
            m_local_to_world_scale: scale,
            m_world_to_local_scale: scale_inv,
            sizes: *local_sizes,
            material: *mat,
            bounding_box: BoundingBox::draw_always(),
        }
    }
}

impl Marched for MarchedTorus {
    fn local_sdf(&self,p: &Point3) -> f32 {
        let p2 = *p;//- center, that is (0,0,0) in local coords
        let q = Vec3::new(Vec3::new(p2.x(),p2.z(),0.).length()-self.sizes.x(),p2.y(),0.);
        return q.length() - self.sizes.y();
    }
    fn material(&self) -> &Material{
        return &self.material;
    }
    fn to_local(&self,p: &Vec4) -> Vec4{
        return self.m_world_to_local_translate_rotate.dot(&(*p*self.m_world_to_local_scale));
    }
    fn to_world(&self,p: &Vec4) -> Vec4{
        return self.m_local_to_world_translate_rotate.dot(p)*self.m_local_to_world_scale;
    }
    fn to_world_f(&self,f: f32) -> f32{//http://jamie-wong.com/2016/07/15/ray-marching-signed-distance-functions/#non-uniform-scaling-and-beyond
        return f*self.m_local_to_world_scale.xyz().min_val();
    }
}
impl Bounded for MarchedTorus {
    fn build_world_bounding_box(&self) -> BoundingBox3D {
        let size = self.sizes.x().max(self.sizes.y());
        let vsize = Vec3::new(size,size,size)*self.m_local_to_world_scale.xyz();
        let bb = BoundingBox3D::new(&-vsize,&vsize).dot(&self.m_local_to_world_translate_rotate);
        return bb;
    }
    fn hit_bounding_box(&self,dir: &Vec3) -> bool{ 
        self.bounding_box.hit(dir) 
    }
    fn get_bounding_box(&self) -> BoundingBox { self.bounding_box }
}