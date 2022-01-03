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
    pub fn to_z1(&self) -> Self {
        *self / self.z()
    }
    pub fn abs(&self) -> Self {
        Vec3::new(self.x().abs(),self.y().abs(),self.z().abs())
    }
    pub fn max_val(&self) -> f32{
        self.x().max(self.y().max(self.z()))
    }
    pub fn max(&self,v: &Vec3) -> Self {
        Vec3::new(self.x().max(v.x()),self.y().max(v.y()),self.z().max(v.z()))
    }
    pub fn min(&self,v: &Vec3) -> Self {
        Vec3::new(self.x().min(v.x()),self.y().min(v.y()),self.z().min(v.z()))
    }
    pub fn sqrt(&self) -> Self{
        Vec3::new(self.x().sqrt(),self.y().sqrt(),self.z().sqrt())
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
    #[inline]
    pub fn to_u8x3(&self) -> (u8,u8,u8){
        return ((self.e[0]*256.0) as u8,(self.e[1]*256.0) as u8,(self.e[2]*256.0) as u8);
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
impl AddAssign<&Vec3> for Vec3{
    fn add_assign(&mut self, other: &Vec3){
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
impl MulAssign<&Vec3> for Vec3{
    fn mul_assign(&mut self, other: &Vec3){
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
    fn mul(self, other: Vec3) -> Vec3{
        return other*self;
    }
}
impl Div<f32> for Vec3{
    type Output = Self;
    fn div(self, scalar: f32) -> Vec3{
        self * (1./scalar)
    }
}