
use crate::math::vec3::{Vec3,Point3};
use crate::math::mat4x4::{Mat4x4};
use crate::utils::INF;

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
    pub fn to_bounding_box(&self,m_to_camera: &Mat4x4) -> BoundingBox{
        let p1 = m_to_camera.dot_p3(&self.minp).to_z1();
        let p2 = m_to_camera.dot_p3(&Point3::new(self.minp.x(),self.minp.y(),self.maxp.z())).to_z1();
        let p3 = m_to_camera.dot_p3(&Point3::new(self.minp.x(),self.maxp.y(),self.minp.z())).to_z1();
        let p4 = m_to_camera.dot_p3(&Point3::new(self.minp.x(),self.maxp.y(),self.maxp.z())).to_z1();
        let p5 = m_to_camera.dot_p3(&Point3::new(self.maxp.x(),self.minp.y(),self.minp.z())).to_z1();
        let p6 = m_to_camera.dot_p3(&Point3::new(self.maxp.x(),self.minp.y(),self.maxp.z())).to_z1();
        let p7 = m_to_camera.dot_p3(&Point3::new(self.maxp.x(),self.maxp.y(),self.minp.z())).to_z1();
        let p8 = m_to_camera.dot_p3(&self.maxp).to_z1();
        let minp = p1.min(&p2.min(&p3.min(&p4.min(&p5.min(&p6.min(&p7.min(&p8)))))));
        let maxp = p1.max(&p2.max(&p3.max(&p4.max(&p5.max(&p6.max(&p7.max(&p8)))))));
        return BoundingBox::new(
            minp.x(),minp.y(),
            maxp.x(),maxp.y(),
        );
    }
}