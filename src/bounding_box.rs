
use crate::math::vec3::{Vec3,Point3};
use crate::math::mat4x4::{Mat4x4};
use crate::utils::INF;
use crate::camera::Camera;

#[derive(Copy,Clone,Debug)]
pub struct BoundingBox {
    pub bottomleft_x: f32,
    pub bottomleft_y: f32,
    pub topright_x: f32,
    pub topright_y: f32,
}

impl BoundingBox {
    pub fn draw_always() -> Self { Self::new(-INF,-INF,INF,INF) }
    pub fn new(bottomleft_x: f32,bottomleft_y: f32,topright_x: f32,topright_y: f32) -> Self { 
        Self{bottomleft_x, bottomleft_y, topright_x, topright_y} 
    }
    pub fn hit(&self,dir: &Vec3) -> bool {
        return self.bottomleft_x <= dir.x() && dir.x() <= self.topright_x 
          &&   self.bottomleft_y <= dir.y() && dir.y() <= self.topright_y; 
    }
}

#[derive(Copy,Clone,Debug)]
pub struct BoundingBox3D {
    pub minp: Point3,
    pub maxp: Point3,
}

impl BoundingBox3D {
    pub fn draw_always() -> Self {Self::new(&Point3::new(-INF,-INF,-INF),&Point3::new(INF,INF,INF))}
    pub fn new(minp: &Point3,maxp: &Point3) -> Self { Self{minp: *minp,maxp: *maxp} }
    pub fn dot(&self,m: &Mat4x4) -> Self {
        let p1 = m.dot_p3(&self.minp);
        let p2 = m.dot_p3(&Point3::new(self.minp.x(),self.minp.y(),self.maxp.z()));
        let p3 = m.dot_p3(&Point3::new(self.minp.x(),self.maxp.y(),self.minp.z()));
        let p4 = m.dot_p3(&Point3::new(self.minp.x(),self.maxp.y(),self.maxp.z()));
        let p5 = m.dot_p3(&Point3::new(self.maxp.x(),self.minp.y(),self.minp.z()));
        let p6 = m.dot_p3(&Point3::new(self.maxp.x(),self.minp.y(),self.maxp.z()));
        let p7 = m.dot_p3(&Point3::new(self.maxp.x(),self.maxp.y(),self.minp.z()));
        let p8 = m.dot_p3(&self.maxp);
        let minp = p1.min(&p2.min(&p3.min(&p4.min(&p5.min(&p6.min(&p7.min(&p8)))))));
        let maxp = p1.max(&p2.max(&p3.max(&p4.max(&p5.max(&p6.max(&p7.max(&p8)))))));
        return BoundingBox3D{minp,maxp};
    }
    pub fn project(&self,cam: &Camera) -> BoundingBox {
        //http://medialab.di.unipi.it/web/IUM/Waterloo/node47.html#SECTION00810000000000000000
        //(x,y,z) -> (xn/z,yn/z,n)
        let p1 = self.minp;
        let p2 = self.maxp;
        let v1 = cam.to_uv(&(p1));
        let v2 = cam.to_uv(&(Point3::new(p1.x(),p1.y(),p2.z())));
        let v3 = cam.to_uv(&(Point3::new(p1.x(),p2.y(),p1.z())));
        let v4 = cam.to_uv(&(Point3::new(p1.x(),p2.y(),p2.z())));
        let v5 = cam.to_uv(&(Point3::new(p2.x(),p1.y(),p1.z())));
        let v6 = cam.to_uv(&(Point3::new(p2.x(),p1.y(),p2.z())));
        let v7 = cam.to_uv(&(Point3::new(p2.x(),p2.y(),p1.z())));
        let v8 = cam.to_uv(&(p2));
        let minu = v1.0.min(v2.0).min(v3.0).min(v4.0).min(v5.0).min(v6.0).min(v7.0).min(v8.0);
        let minv = v1.1.min(v2.1).min(v3.1).min(v4.1).min(v5.1).min(v6.1).min(v7.1).min(v8.1);
        let maxu = v1.0.max(v2.0).max(v3.0).max(v4.0).max(v5.0).max(v6.0).max(v7.0).max(v8.0);
        let maxv = v1.1.max(v2.1).max(v3.1).max(v4.1).max(v5.1).max(v6.1).max(v7.1).max(v8.1);
        return BoundingBox::new(minu,minv,maxu,maxv);
    }
}

pub trait Bounded {
    fn build_world_bounding_box(&self) -> BoundingBox3D { BoundingBox3D::draw_always() }
    fn get_bounding_box(&self) -> BoundingBox { BoundingBox::draw_always() }
    fn hit_bounding_box(&self,_dir: &Vec3) -> bool{ true }
}