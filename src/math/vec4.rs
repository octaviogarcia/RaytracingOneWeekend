#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct Vec4 {
    pub e: [f32;4]
}

impl std::fmt::Display for Vec4 {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "Vec4{{x: {}, y: {}, z: {},w: {}}}", self.x(), self.y(), self.z(), self.w())
    }
}

use crate::utils::{max,min};
use super::vec3::Vec3;

impl Vec4{
    pub const ZERO : Self = Self{ e: [0.,0.,0.,0.]};
    pub fn new(x: f32,y: f32,z: f32,w: f32) -> Self { Self{e: [x,y,z,w]} }
    pub fn new_v3(v: &Vec3,w: f32) -> Self { Self {e: [v.x(),v.y(),v.z(),w]} }
    pub fn new_unit(x: f32,y: f32,z: f32,w: f32) -> Self { Self{e: [x,y,z,w]}.unit() }

    pub fn x(&self) -> f32{ self.e[0] }
    pub fn y(&self) -> f32{ self.e[1] }
    pub fn z(&self) -> f32{ self.e[2] }
    pub fn w(&self) -> f32{ self.e[3] }

    pub fn dot(&self,other: Self) -> f32{
        self.e[0]*other.e[0]+self.e[1]*other.e[1]+self.e[2]*other.e[2]+self.e[3]*other.e[3]
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
        Vec4::new(self.x().abs(),self.y().abs(),self.z().abs(),self.w().abs())
    }
    pub fn max(&self,v: &Vec4) -> Self {
        Vec4::new(max(self.x(),v.x()),max(self.y(),v.y()),max(self.z(),v.z()),max(self.w(),v.w()))
    }
    pub fn min(&self,v: &Vec4) -> Self {
        Vec4::new(min(self.x(),v.x()),min(self.y(),v.y()),min(self.z(),v.z()),min(self.w(),v.w()))
    }
}
use std::ops::{Neg,Index,IndexMut,AddAssign,SubAssign,MulAssign,DivAssign,Add,Sub,Mul,Div};

impl Neg for Vec4{
    type Output = Self;
    fn neg(self) -> Vec4{ Vec4{ e : [-self.e[0],-self.e[1],-self.e[2],-self.e[3]]} }
}
impl Index<usize> for Vec4{
    type Output = f32;
    fn index(&self, idx: usize) -> &Self::Output  { &self.e[idx]  }
}
impl IndexMut<usize> for Vec4{
    fn index_mut(&mut self, idx: usize) -> &mut Self::Output { &mut self.e[idx]  }
}
impl AddAssign<Vec4> for Vec4{
    fn add_assign(&mut self, other: Vec4){
        *self = Vec4 { e: [
            self.e[0]+other.e[0],
            self.e[1]+other.e[1],
            self.e[2]+other.e[2],
            self.e[3]+other.e[3],
        ]}
    }
}
impl SubAssign<Vec4> for Vec4{
    fn sub_assign(&mut self, other: Vec4){
        *self += -other;
    }
}
impl MulAssign<Vec4> for Vec4{
    fn mul_assign(&mut self, other: Vec4){
        *self = Vec4 { e: [
            self.e[0]*other.e[0],
            self.e[1]*other.e[1],
            self.e[2]*other.e[2],
            self.e[3]*other.e[3],
        ]}
    }
}
impl MulAssign<f32> for Vec4{
    fn mul_assign(&mut self, scalar: f32){
        *self = Vec4 { e: [
            self.e[0]*scalar,
            self.e[1]*scalar,
            self.e[2]*scalar,
            self.e[3]*scalar,
        ]}
    }
}
impl DivAssign<f32> for Vec4{
    fn div_assign(&mut self, scalar: f32){
        *self *= 1./scalar;
    }
}
impl Add<Vec4> for Vec4{
    type Output = Self;
    fn add(self, other: Vec4) -> Vec4{
        Vec4 { e : [
            self.e[0]+other.e[0],
            self.e[1]+other.e[1],
            self.e[2]+other.e[2],
            self.e[3]+other.e[3],
        ]}
    }
}
impl Sub<Vec4> for Vec4{
    type Output = Self;
    fn sub(self, other: Vec4) -> Vec4{
        self + (-other)
    }
}
impl Mul<Vec4> for Vec4{
    type Output = Self;
    fn mul(self, other: Vec4) -> Vec4{
        Vec4 { e : [
            self.e[0]*other.e[0],
            self.e[1]*other.e[1],
            self.e[2]*other.e[2],
            self.e[3]*other.e[3],
        ]}
    }
}
impl Div<Vec4> for Vec4{
    type Output = Self;
    fn div(self, other: Vec4) -> Vec4{
        Vec4 { e : [
            self.e[0]/other.e[0],
            self.e[1]/other.e[1],
            self.e[2]/other.e[2],
            self.e[3]/other.e[3],
        ]}
    }
}
impl Mul<f32> for Vec4{
    type Output = Self;
    fn mul(self, scalar: f32) -> Vec4{
        Vec4 { e : [
            self.e[0]*scalar,
            self.e[1]*scalar,
            self.e[2]*scalar,
            self.e[3]*scalar,
        ]}
    }
}
impl Mul<Vec4> for f32{
    type Output = Vec4;
    fn mul(self, other: Vec4) -> Vec4{
        return other*self;
    }
}
impl Div<f32> for Vec4{
    type Output = Self;
    fn div(self, scalar: f32) -> Vec4{
        self * (1./scalar)
    }
}