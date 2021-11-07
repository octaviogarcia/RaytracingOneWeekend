use super::vec3::Vec3;
//(-e[0]-)
//(-e[1]-)
//(-e[2]-)
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct Mat3x3 {
    pub e: [Vec3;3]
}

impl Mat3x3 {
    #[allow(dead_code)]
    pub const ZERO : Self     = Self{ e: [Vec3{e:[0.,0.,0.]},Vec3{e:[0.,0.,0.]},Vec3{e:[0.,0.,0.]}]};
    #[allow(dead_code)]
    pub const IDENTITY : Self = Self{ e: [Vec3{e:[1.,0.,0.]},Vec3{e:[0.,1.,0.]},Vec3{e:[0.,0.,1.]}]};
    #[allow(dead_code)]
    pub fn new_9f(e00: f32,e01: f32,e02: f32,e10: f32,e11: f32,e12: f32,e20: f32,e21: f32,e22: f32) -> Self{
        return Self{e: [Vec3::new(e00,e01,e02),Vec3::new(e10,e11,e12),Vec3::new(e20,e21,e22)]};
    }
    pub fn new_3vec(v0: &Vec3,v1: &Vec3,v2: &Vec3) -> Self{
        return Self{e: [*v0,*v1,*v2]};
    }
    pub fn new_3vec_vert(v0: &Vec3,v1: &Vec3,v2: &Vec3) -> Self{
        return Self{e: [Vec3::new(v0.x(),v1.y(),v2.x()),
                        Vec3::new(v0.y(),v1.y(),v2.y()),
                        Vec3::new(v0.z(),v1.z(),v2.z())]};
    }
    pub fn dot(&self,v: &Vec3) -> Vec3{
        return Vec3::new(self.e[0].dot(*v),self.e[1].dot(*v),self.e[2].dot(*v));
    }
    #[allow(dead_code)]
    pub fn dot_mat(&self,m: &Self) -> Self{
        let t = m.transpose();
        return Self::new_9f(self.at_row(0).dot(t.at_row(0)),self.at_row(0).dot(t.at_row(1)),self.at_row(0).dot(t.at_row(2)),
                           self.at_row(1).dot(t.at_row(0)),self.at_row(1).dot(t.at_row(1)),self.at_row(1).dot(t.at_row(2)),
                           self.at_row(2).dot(t.at_row(0)),self.at_row(2).dot(t.at_row(1)),self.at_row(2).dot(t.at_row(2)));
    }
    pub fn at(&self,row: usize,col: usize) -> &f32{
        return &self.e[row].e[col];
    }
    #[allow(dead_code)]
    pub fn at_row(&self,row: usize) -> Vec3{
        return self.e[row];
    }
    #[allow(dead_code)]
    pub fn at_col(&self,col: usize) -> Vec3{
        return Vec3::new(self.e[0][col],self.e[1][col],self.e[2][col]);
    }
    #[allow(dead_code)]
    pub fn transpose(&self) -> Mat3x3{
        return Self::new_3vec(&self.at_col(0),&self.at_col(1),&self.at_col(2));   
    }
    pub fn determinant(&self) -> f32{
        let mut sum: f32 = 0.;//Kahan summation
        let mut c: f32 = 0.;
        let input = [(self.at(0,0)*self.at(1,1)*self.at(2,2)),
                     (self.at(0,1)*self.at(1,2)*self.at(2,0)),
                     (self.at(0,2)*self.at(1,0)*self.at(2,1)),
                     (-self.at(0,0)*self.at(1,2)*self.at(2,1)),
                     (-self.at(0,1)*self.at(1,0)*self.at(2,2)),
                     (-self.at(0,2)*self.at(1,1)*self.at(2,0))];
        for i in input{
            let y = i - c;
            let t = sum + y;
            c = (t - sum) - y;
            sum = t;
        }
        return sum;
    }
    pub fn inverse(&self) -> Mat3x3{
        let det = self.determinant();
        let row0 = Vec3::new(self.at(1,1)*self.at(2,2)-self.at(1,2)*self.at(2,1),
                             self.at(0,2)*self.at(2,1)-self.at(0,1)*self.at(2,2),
                             self.at(0,1)*self.at(1,2)-self.at(0,2)*self.at(1,1));
        let row1 = Vec3::new(self.at(1,2)*self.at(2,0)-self.at(1,0)*self.at(2,2),
                             self.at(0,0)*self.at(2,2)-self.at(0,2)*self.at(2,0),
                             self.at(0,2)*self.at(1,0)-self.at(0,0)*self.at(1,2));
        let row2 = Vec3::new(self.at(1,0)*self.at(2,1)-self.at(1,1)*self.at(2,0),
                             self.at(0,1)*self.at(2,0)-self.at(0,0)*self.at(2,1),
                             self.at(0,0)*self.at(1,1)-self.at(0,1)*self.at(1,0));
        return Self::new_3vec(&row0,&row1,&row2)/det;
    }
}

use std::ops::{Mul,Div};

impl Mul<f32> for Mat3x3{
    type Output = Self;
    fn mul(self, scalar: f32) -> Self{
        return Self::new_3vec(&(self.e[0]*scalar),&(self.e[1]*scalar),&(self.e[2]*scalar));
    }
}
impl Div<f32> for Mat3x3{
    type Output = Self;
    fn div(self, scalar: f32) -> Self{
        self * (1./scalar)
    }
}