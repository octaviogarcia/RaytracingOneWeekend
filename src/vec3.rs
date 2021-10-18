#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct Vec3 {
    pub e: [f32;3]
}

impl std::fmt::Display for Vec3 {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "Vec3{{x: {}, y: {}, z: {}}}", self.x(), self.y(), self.z())
    }
}

// Type aliases for vec3
pub type Point3 = Vec3;   // 3D point
pub type Color  = Vec3;   // RGB color
pub type UnitVec3 = Vec3; //Type alias just to document unit vectors

use crate::utils;
use utils::MyRandom;
use utils::{max,min,abs};

impl Vec3{
    pub const ZERO : Self = Self{ e: [0.,0.,0.]};
    pub fn new(x: f32,y: f32,z: f32) -> Self { Self{e: [x,y,z]} }
    pub fn new_unit(x: f32,y: f32,z: f32) -> Self { Self{e: [x,y,z]}.unit() }

    pub fn x(&self) -> f32{ self.e[0] }
    pub fn y(&self) -> f32{ self.e[1] }
    pub fn z(&self) -> f32{ self.e[2] }

    pub fn dot(&self,other: Self) -> f32{
        self.e[0]*other.e[0]+self.e[1]*other.e[1]+self.e[2]*other.e[2]
    }
    pub fn length_squared(&self) -> f32 {
        self.dot(*self)
    }
    pub fn length(&self) -> f32 {
        self.length_squared().sqrt()
    }
    pub fn unit(&self) -> Self {
        *self / self.length()
    }
    pub fn abs(&self) -> Self {
        Vec3::new(abs(self.x()),abs(self.y()),abs(self.z()))
    }
    pub fn max(&self,v: &Vec3) -> Self {
        Vec3::new(max(self.x(),v.x()),max(self.y(),v.y()),max(self.z(),v.z()))
    }
    pub fn min(&self,v: &Vec3) -> Self {
        Vec3::new(min(self.x(),v.x()),min(self.y(),v.y()),min(self.z(),v.z()))
    }
    pub fn cross(&self, other: Self) -> Self{
        Self{ e : [
            self.e[1]*other.e[2] - self.e[2]*other.e[1],
            self.e[2]*other.e[0] - self.e[0]*other.e[2],
            self.e[0]*other.e[1] - self.e[1]*other.e[0],
        ]}
    }
    pub fn near_zero(&self) -> bool {
        let eps = 1e-8;
        return (self.e[0].abs() < eps) && (self.e[1].abs() < eps) && (self.e[2].abs() < eps);
    }
}

//Random implementations
impl MyRandom for Vec3{
    fn rand() -> Self { Self{e: [f32::rand(),f32::rand(),f32::rand()]} }
    fn rand_range(fmin: f32,fmax: f32) -> Self {
        Self{e: [
            f32::rand_range(fmin,fmax),
            f32::rand_range(fmin,fmax),
            f32::rand_range(fmin,fmax),
        ]}
    }
}

impl Vec3{
    pub fn rand_in_unit_sphere() -> Self{
        loop {
            let p = Self::rand_range(-1.,1.);
            if p.length_squared() < 1. {return p;};
        }
    }
    pub fn rand_unit_vector() -> Self{
        return Self::rand_in_unit_sphere().unit();
    }
    pub fn rand_in_hemisphere(normal: &Self) -> Self{
        let in_unit_sphere = Self::rand_in_unit_sphere();
        if in_unit_sphere.dot(*normal) > 0. {//Same hemisphere
            return in_unit_sphere;
        }
        return -in_unit_sphere;//Invert it
    }
    pub fn rand_in_unit_disc() -> Self{
        loop {
            let p = Self::new(f32::rand_range(-1.,1.),f32::rand_range(-1.,1.),0.);
            if p.length_squared() < 1. {return p;};
        }
    }
}

use std::ops::{Neg,Index,IndexMut,AddAssign,SubAssign,MulAssign,DivAssign,Add,Sub,Mul,Div};

impl Neg for Vec3{
    type Output = Self;
    fn neg(self) -> Vec3{ Vec3{ e : [-self.e[0],-self.e[1],-self.e[2]]} }
}
impl Index<usize> for Vec3{
    type Output = f32;
    fn index(&self, idx: usize) -> &Self::Output  { &self.e[idx]  }
}
impl IndexMut<usize> for Vec3{
    fn index_mut(&mut self, idx: usize) -> &mut Self::Output { &mut self.e[idx]  }
}
impl AddAssign<Vec3> for Vec3{
    fn add_assign(&mut self, other: Vec3){
        *self = Vec3 { e: [
            self.e[0]+other.e[0],
            self.e[1]+other.e[1],
            self.e[2]+other.e[2],
        ]}
    }
}
impl SubAssign<Vec3> for Vec3{
    fn sub_assign(&mut self, other: Vec3){
        *self += -other;
    }
}
impl MulAssign<Vec3> for Vec3{
    fn mul_assign(&mut self, other: Vec3){
        *self = Vec3 { e: [
            self.e[0]*other.e[0],
            self.e[1]*other.e[1],
            self.e[2]*other.e[2],
        ]}
    }
}
impl MulAssign<f32> for Vec3{
    fn mul_assign(&mut self, scalar: f32){
        *self = Vec3 { e: [
            self.e[0]*scalar,
            self.e[1]*scalar,
            self.e[2]*scalar,
        ]}
    }
}
impl DivAssign<f32> for Vec3{
    fn div_assign(&mut self, scalar: f32){
        *self *= 1./scalar;
    }
}
impl Add<Vec3> for Vec3{
    type Output = Self;
    fn add(self, other: Vec3) -> Vec3{
        Vec3 { e : [
            self.e[0]+other.e[0],
            self.e[1]+other.e[1],
            self.e[2]+other.e[2],
        ]}
    }
}
impl Sub<Vec3> for Vec3{
    type Output = Self;
    fn sub(self, other: Vec3) -> Vec3{
        self + (-other)
    }
}
impl Mul<Vec3> for Vec3{
    type Output = Self;
    fn mul(self, other: Vec3) -> Vec3{
        Vec3 { e : [
            self.e[0]*other.e[0],
            self.e[1]*other.e[1],
            self.e[2]*other.e[2],
        ]}
    }
}
impl Div<Vec3> for Vec3{
    type Output = Self;
    fn div(self, other: Vec3) -> Vec3{
        Vec3 { e : [
            self.e[0]/other.e[0],
            self.e[1]/other.e[1],
            self.e[2]/other.e[2],
        ]}
    }
}
impl Mul<f32> for Vec3{
    type Output = Self;
    fn mul(self, scalar: f32) -> Vec3{
        Vec3 { e : [
            self.e[0]*scalar,
            self.e[1]*scalar,
            self.e[2]*scalar,
        ]}
    }
}
impl Mul<Vec3> for f32{
    type Output = Vec3;
    fn mul(self, vec3: Vec3) -> Vec3{
        return vec3*self;
    }
}
impl Div<f32> for Vec3{
    type Output = Self;
    fn div(self, scalar: f32) -> Vec3{
        self * (1./scalar)
    }
}

//(-e[0]-)
//(-e[1]-)
//(-e[2]-)
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct Mat3x3 {
    pub e: [Vec3;3]
}

impl Mat3x3 {
    pub const ZERO : Self     = Self{ e: [Vec3{e:[0.,0.,0.]},Vec3{e:[0.,0.,0.]},Vec3{e:[0.,0.,0.]}]};
    pub const IDENTITY : Self = Self{ e: [Vec3{e:[1.,0.,0.]},Vec3{e:[0.,1.,0.]},Vec3{e:[0.,0.,1.]}]};
    pub fn new9f(e00: f32,e01: f32,e02: f32,e10: f32,e11: f32,e12: f32,e20: f32,e21: f32,e22: f32) -> Self{
        return Self{e: [Vec3::new(e00,e01,e02),Vec3::new(e10,e11,e12),Vec3::new(e20,e21,e22)]};
    }
    pub fn new3v(v0: &Vec3,v1: &Vec3,v2: &Vec3) -> Self{
        return Self{e: [*v0,*v1,*v2]};
    }
    pub fn new3v_vert(v0: &Vec3,v1: &Vec3,v2: &Vec3) -> Self{
        return Self{e: [Vec3::new(v0.x(),v1.y(),v2.x()),
                        Vec3::new(v0.y(),v1.y(),v2.y()),
                        Vec3::new(v0.z(),v1.z(),v2.z())]};
    }
    pub fn dot(&self,v: &Vec3) -> Vec3{
        return Vec3::new(self.e[0].dot(*v),self.e[1].dot(*v),self.e[2].dot(*v));
    }
    pub fn dot_mat(&self,v: &Self) -> Self{
        let t = self.transpose();
        return Self::new9f(self.at_row(0).dot(t.at_row(0)),self.at_row(0).dot(t.at_row(1)),self.at_row(0).dot(t.at_row(2)),
                           self.at_row(1).dot(t.at_row(0)),self.at_row(1).dot(t.at_row(1)),self.at_row(1).dot(t.at_row(2)),
                           self.at_row(2).dot(t.at_row(0)),self.at_row(2).dot(t.at_row(1)),self.at_row(2).dot(t.at_row(2)));
    }
    pub fn at(&self,row: usize,col: usize) -> f32{
        return self.e[row].e[col];
    }
    pub fn at_row(&self,row: usize) -> Vec3{
        return self.e[row];
    }
    pub fn at_col(&self,col: usize) -> Vec3{
        return Vec3::new(self.e[0][col],self.e[1][col],self.e[2][col]);
    }
    pub fn transpose(&self) -> Mat3x3{
        return Self::new3v(&self.at_col(0),&self.at_col(1),&self.at_col(2));   
    }
    pub fn determinant(&self) -> f32{
        return self.at(0,0)*self.at(1,1)*self.at(2,2)
             + self.at(0,1)*self.at(1,2)*self.at(2,0)
             + self.at(0,2)*self.at(1,0)*self.at(2,1)
             - self.at(0,0)*self.at(1,2)*self.at(2,1)
             - self.at(0,1)*self.at(1,0)*self.at(2,2)
             - self.at(0,2)*self.at(1,1)*self.at(2,0);
    }
    pub fn inverse(&self) -> Mat3x3{
        let row0 = Vec3::new(self.at(1,1)*self.at(2,2)-self.at(1,2)*self.at(2,1),
                             self.at(0,2)*self.at(2,1)-self.at(0,1)*self.at(2,2),
                             self.at(0,1)*self.at(1,2)-self.at(0,2)*self.at(1,1));
        let row1 = Vec3::new(self.at(1,2)*self.at(2,0)-self.at(1,0)*self.at(2,2),
                             self.at(0,0)*self.at(2,2)-self.at(0,2)*self.at(2,0),
                             self.at(0,2)*self.at(1,0)-self.at(0,0)*self.at(1,2));
        let row2 = Vec3::new(self.at(1,0)*self.at(2,1)-self.at(1,1)*self.at(2,0),
                             self.at(0,1)*self.at(2,0)-self.at(0,0)*self.at(2,1),
                             self.at(0,0)*self.at(1,1)-self.at(0,1)*self.at(1,0));
        return Self::new3v(&row0,&row1,&row2)/self.determinant();
    }
}
impl Mul<f32> for Mat3x3{
    type Output = Self;
    fn mul(self, scalar: f32) -> Self{
        return Self::new3v(&(self.e[0]*scalar),&(self.e[1]*scalar),&(self.e[2]*scalar));
    }
}
impl Div<f32> for Mat3x3{
    type Output = Self;
    fn div(self, scalar: f32) -> Self{
        self * (1./scalar)
    }
}