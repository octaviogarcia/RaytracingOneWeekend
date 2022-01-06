use super::vec4::Vec4;
use super::vec3::{Vec3,Point3};
use super::mat3x3::Mat3x3;
//(-e[0]-)
//(-e[1]-)
//(-e[2]-)
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct Mat4x4 {
    pub e: [Vec4;4]
}

impl Mat4x4 {
    #[allow(dead_code)]
    pub const ZERO : Self     = Self{ e: [Vec4::ZERO,Vec4::ZERO,Vec4::ZERO,Vec4::ZERO]};
    #[allow(dead_code)]
    pub const IDENTITY : Self = Self{ e: [Vec4{e:[1.,0.,0.,0.]},Vec4{e:[0.,1.,0.,0.]},Vec4{e:[0.,0.,1.,0.]},Vec4{e:[0.,0.,0.,1.]}]};
    #[allow(dead_code)]
    pub fn new_16f(e00: f32,e01: f32,e02: f32,e03: f32,
                  e10: f32,e11: f32,e12: f32,e13: f32,
                  e20: f32,e21: f32,e22: f32,e23: f32,
                  e30: f32,e31: f32,e32: f32,e33: f32) -> Self{
        Self{e: [Vec4::new(e00,e01,e02,e03),
                 Vec4::new(e10,e11,e12,e13),
                 Vec4::new(e20,e21,e22,e23),
                 Vec4::new(e30,e31,e32,e33)]}
    }
    #[allow(dead_code)]
    pub fn new_4vec(v0: &Vec4,v1: &Vec4,v2: &Vec4,v3: &Vec4) -> Self{
        Self{e: [*v0,*v1,*v2,*v3]}
    }
    #[allow(dead_code)]
    pub fn new_4vec_vert(v0: &Vec4,v1: &Vec4,v2: &Vec4,v3: &Vec4) -> Self{
        Self{e: [Vec4::new(v0.x(),v1.x(),v2.x(),v3.x()),
                 Vec4::new(v0.y(),v1.y(),v2.y(),v3.y()),
                 Vec4::new(v0.z(),v1.z(),v2.z(),v3.z()),
                 Vec4::new(v0.w(),v1.w(),v2.w(),v3.w())]}
    }
    pub fn new_m3x3(m: &Mat3x3) -> Self{
        Self{e: [Vec4::new_v3(&m.e[0]),
                 Vec4::new_v3(&m.e[1]),
                 Vec4::new_v3(&m.e[2]),
                 Vec4::new(0.,0.,0.,1.)]}
    }
    #[allow(dead_code)]
    pub fn dot(&self,v: &Vec4) -> Vec4{
        Vec4::new(self.e[0].dot(*v),self.e[1].dot(*v),self.e[2].dot(*v),self.e[3].dot(*v))
    }
    #[allow(dead_code)]
    pub fn dot_p3(&self,v: &Point3) -> Point3{
        let v4 = Vec4::new_p3(v);
        Point3::new(self.e[0].dot(v4),self.e[1].dot(v4),self.e[2].dot(v4))
    }
    #[allow(dead_code)]
    pub fn dot_v3(&self,v: &Vec3) -> Vec3{
        let v4 = Vec4::new_v3(v);
        Vec3::new(self.e[0].dot(v4),self.e[1].dot(v4),self.e[2].dot(v4))
    }
    #[allow(dead_code)]
    pub fn fast_homogenous_inverse(&self) -> Self{//https://stackoverflow.com/questions/155670/invert-4x4-matrix-numerical-most-stable-solution-needed
        //M = TS -> M^-1 = (TS)^-1 = S^-1 T^-1
        let s_inv = Self::new_m3x3(&Mat3x3::new_3vec(&self.e[0].xyz(),&self.e[1].xyz(),&self.e[2].xyz()).inverse());
        let t_inv = Self::new_translate(&-self.at_col(3).xyz());
        return s_inv.dot_mat(&t_inv);
    }
    #[allow(dead_code)]
    pub fn dot_mat(&self,m: &Self) -> Self{
        let t = m.transpose();
        return Self::new_16f(self.at_row(0).dot(t.at_row(0)),self.at_row(0).dot(t.at_row(1)),self.at_row(0).dot(t.at_row(2)),self.at_row(0).dot(t.at_row(3)),
                            self.at_row(1).dot(t.at_row(0)),self.at_row(1).dot(t.at_row(1)),self.at_row(1).dot(t.at_row(2)),self.at_row(1).dot(t.at_row(3)),
                            self.at_row(2).dot(t.at_row(0)),self.at_row(2).dot(t.at_row(1)),self.at_row(2).dot(t.at_row(2)),self.at_row(2).dot(t.at_row(3)),
                            self.at_row(3).dot(t.at_row(0)),self.at_row(3).dot(t.at_row(1)),self.at_row(3).dot(t.at_row(2)),self.at_row(3).dot(t.at_row(3)));
    }
    #[allow(dead_code)]      
    pub fn at(&self,row: usize,col: usize) -> f32{
        self.e[row].e[col]
    }
    #[allow(dead_code)]
    pub fn at_row(&self,row: usize) -> Vec4{
        self.e[row]
    }
    pub fn at_col(&self,col: usize) -> Vec4{
        Vec4::new(self.e[0][col],self.e[1][col],self.e[2][col],self.e[3][col])
    }
    #[allow(dead_code)]
    pub fn transpose(&self) -> Self{
        Self::new_4vec(&self.at_col(0),&self.at_col(1),&self.at_col(2),&self.at_col(3))
    }
    pub fn new_translate(v: &Vec3) -> Self{
        Self{ e: [Vec4{e:[1.,0.,0.,v.x()]},
                  Vec4{e:[0.,1.,0.,v.y()]},
                  Vec4{e:[0.,0.,1.,v.z()]},
                  Vec4{e:[0.,0.,0.,   1.]}] }
    }
    pub fn new_scale(v: &Vec3) -> Self{
        Self{ e: [ Vec4{e:[v.x(),   0.,    0.,0.]},
                   Vec4{e:[   0.,v.y(),    0.,0.]},
                   Vec4{e:[   0.,   0., v.z(),0.]},
                   Vec4{e:[   0.,   0.,    0.,1.]}] }
    }
    pub fn new_rotate_x(f: f32) -> Self{
        let c = f.cos();
        let s = f.sin();
        Self{ e: [Vec4{e:[1.,0.,0.,0.]},
                  Vec4{e:[0., c, s,0.]},
                  Vec4{e:[0.,-s, c,0.]},
                  Vec4{e:[0.,0.,0.,1.]}] }
    }
    pub fn new_rotate_y(f: f32) -> Self{
        let c = f.cos();
        let s = f.sin();
        Self{ e: [Vec4{e:[ c,0.,-s,0.]},
                  Vec4{e:[0.,1.,0.,0.]},
                  Vec4{e:[ s,0., c,0.]},
                  Vec4{e:[0.,0.,0.,1.]}] }
    }
    pub fn new_rotate_z(f: f32) -> Self{
        let c = f.cos();
        let s = f.sin();
        Self{ e: [Vec4{e:[ c,-s,0.,0.]},
                  Vec4{e:[ s, c,0.,0.]},
                  Vec4{e:[0.,0.,1.,0.]},
                  Vec4{e:[0.,0.,0.,1.]}] }
    }
    //Decomposes into TRS, so if you invert them individually you need to apply like this S^-1 R^-1 T^-1
    pub fn decompose_into_translate_rotate_scale(&self) -> (Self,Self,Self){
        //@WARNING: Doesn't work if scaling is negative (or shearing... I think)!!! maybe assert() on the auxiliary creation functions
        let a = self.at_col(0).xyz();
        let b = self.at_col(1).xyz();
        let c = self.at_col(2).xyz();
        let d = self.at_col(3).xyz();
        let al = a.length();
        let bl = b.length();
        let cl = c.length();
        let translate = Mat4x4::new_4vec_vert(
            &Vec4::new(1.,0.,0.,0.),&Vec4::new(0.,1.,0.,0.),&Vec4::new(0.,0.,1.,0.),&Vec4::new_p3(&d)
        );
        let rotate    = Mat4x4::new_4vec_vert(
            &Vec4::new_v3(&(a/al)),&Vec4::new_v3(&(b/bl)),&Vec4::new_v3(&(c/cl)),&Vec4::new_p3(&Vec3::ZERO)
        );
        let scale     = Mat4x4::new_4vec(
            &Vec4::new(al,0.,0.,0.),
            &Vec4::new(0.,bl,0.,0.),
            &Vec4::new(0.,0.,cl,0.),
            &Vec4::new(0.,0.,0.,1.),
        ); 
        return (translate,rotate,scale);
    }
    pub fn diag(&self) -> Vec4{
        return Vec4::new(self.at_row(0).x(),self.at_row(1).y(),self.at_row(2).z(),self.at_row(3).w());
    }
}

use std::ops::{Mul,Div,BitXor};

impl Mul<f32> for Mat4x4{
    type Output = Self;
    fn mul(self, scalar: f32) -> Self{
        return Self::new_4vec(&(self.e[0]*scalar),&(self.e[1]*scalar),&(self.e[2]*scalar),&(self.e[3]*scalar));
    }
}
impl Div<f32> for Mat4x4{
    type Output = Self;
    fn div(self, scalar: f32) -> Self{
        self * (1./scalar)
    }
}
impl BitXor<&Mat4x4> for Mat4x4{
    type Output = Self;
    fn bitxor(self,other: &Self) -> Self{
        self.dot_mat(other)
    }
}
impl BitXor<Mat4x4> for Mat4x4{
    type Output = Self;
    fn bitxor(self,other: Self) -> Self{
        self.dot_mat(&other)
    }
}

#[macro_export]
macro_rules! m4x4 {
    (RX $f:expr) => {
        Mat4x4::new_rotate_x($f)
    };
    (RY $f:expr) => {
        Mat4x4::new_rotate_y($f)
    };
    (RZ $f:expr) => {
        Mat4x4::new_rotate_z($f)
    };
    (TR $f1:expr,$f2:expr,$f3:expr) => {
        Mat4x4::new_translate(&Vec3::new($f1,$f2,$f3))
    };
    (TR $v:expr) => {
        Mat4x4::new_translate(&$v)
    };
    (SC $f1:expr,$f2:expr,$f3:expr) => {
        Mat4x4::new_scale(&Vec3::new($f1,$f2,$f3))
    };
    (SC $v:expr) => {
        Mat4x4::new_scale(&$v)
    };
    (ID) => {
        Mat4x4::IDENTITY
    };
    () => {
        Mat4x4::IDENTITY
    };
}
/*
#[macro_export]
macro_rules! transforms {
    ($($e: expr;)*) => {
        m4x4!($e)
    };
}*/
/*
#[macro_export]
macro_rules! transforms {
    ($($ty:ident ; $f:expr);*) => {
        Mat4x4::IDENTITY(^m4x4!($ty $f))*
    };
}
macro_rules! transform {
    ($($traced_ident:ident ; $traced:ty),* | $($marched_ident:ident ; $marched:ty),*) => {
    }
}*/