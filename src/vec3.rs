#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct Vec3 {
    pub e: [f64;3]
}
// Type aliases for vec3
pub type Point3 = Vec3;   // 3D point
pub type Color  = Vec3;   // RGB color

impl Vec3{
    pub const ZERO : Self = Self{ e: [0.,0.,0.]};
    pub fn new(x: f64,y: f64,z: f64) -> Self { Self{e: [x,y,z]} }

    pub fn x(&self) -> &f64{ &self.e[0] }
    pub fn y(&self) -> &f64{ &self.e[1] }
    pub fn z(&self) -> &f64{ &self.e[2] }
    pub fn dot(&self,other: Self) -> f64{
        self.e[0]*other.e[0]+self.e[1]*other.e[1]+self.e[2]*other.e[2]
    }
    pub fn length_squared(&self) -> f64 {
        self.dot(*self)
    }
    pub fn length(&self) -> f64 {
        self.length_squared().sqrt()
    }
    pub fn unit(&self) -> Self {
        *self / self.length()
    }
    pub fn cross(&self, other: Self) -> Self{
        Self{ e : [
            self.e[1]*other.e[2] - self.e[2]*other.e[1],
            self.e[2]*other.e[0] - self.e[0]*other.e[2],
            self.e[0]*other.e[1] - self.e[1]*other.e[0],
        ]}
    }
}

use std::ops::{Neg,Index,IndexMut,AddAssign,SubAssign,MulAssign,DivAssign,Add,Sub,Mul,Div};

impl Neg for Vec3{
    type Output = Self;
    fn neg(self) -> Vec3{ Vec3{ e : [-self.e[0],-self.e[1],-self.e[2]]} }
}
impl Index<usize> for Vec3{
    type Output = f64;
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
impl MulAssign<f64> for Vec3{
    fn mul_assign(&mut self, scalar: f64){
        *self = Vec3 { e: [
            self.e[0]*scalar,
            self.e[1]*scalar,
            self.e[2]*scalar,
        ]}
    }
}
impl DivAssign<f64> for Vec3{
    fn div_assign(&mut self, scalar: f64){
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
impl Mul<f64> for Vec3{
    type Output = Self;
    fn mul(self, scalar: f64) -> Vec3{
        Vec3 { e : [
            self.e[0]*scalar,
            self.e[1]*scalar,
            self.e[2]*scalar,
        ]}
    }
}
impl Mul<Vec3> for f64{
    type Output = Vec3;
    fn mul(self, vec3: Vec3) -> Vec3{
        return vec3*self;
    }
}
impl Div<f64> for Vec3{
    type Output = Self;
    fn div(self, scalar: f64) -> Vec3{
        self * (1./scalar)
    }
}