
use crate::utils::{INF,clamp};
use crate::bounding_box::BoundingBox;

pub trait CellEmptyInitializable: Copy {
    fn empty(&mut self) -> ();
}

const CAMERA_HASH_SIZE: usize = 100;
#[derive(Debug)]
pub struct CameraHash<Cell>{
    pub cells: [[Cell;CAMERA_HASH_SIZE];CAMERA_HASH_SIZE],
    pub bottomleft_x: f32,
    pub bottomleft_y: f32,
    pub topright_x: f32,
    pub topright_y: f32,
}

impl <Cell: CellEmptyInitializable> CameraHash<Cell>{
    pub fn empty(&mut self){
        for i in 0..CAMERA_HASH_SIZE{
            for j in 0..CAMERA_HASH_SIZE{
                self.cells[i][j].empty();
            }
        }
        self.bottomleft_x =  INF;
        self.bottomleft_y =  INF;
        self.topright_x   = -INF;
        self.topright_y   = -INF;
    }
    pub fn set_borders(&mut self,bb: &BoundingBox) -> (){
        if bb.bottomleft_x.is_finite(){
            self.bottomleft_x = bb.bottomleft_x.min(self.bottomleft_x);
        }
        if bb.bottomleft_y.is_finite(){
            self.bottomleft_y = bb.bottomleft_y.min(self.bottomleft_y);
        }
        if bb.topright_x.is_finite(){
            self.topright_x   = bb.topright_x.max(self.topright_x);
        }
        if bb.topright_y.is_finite(){
            self.topright_y   = bb.topright_y.max(self.topright_y);
        }
    }
    pub fn get_indexes(&self,bb: &BoundingBox) -> (usize,usize,usize,usize){
        let left  = (bb.bottomleft_x - self.bottomleft_x)/(self.topright_x-self.bottomleft_x);
        let right = (bb.topright_x   - self.bottomleft_x)/(self.topright_x-self.bottomleft_x);
        let down  = (bb.bottomleft_y - self.bottomleft_y)/(self.topright_y-self.bottomleft_y);
        let up    = (bb.topright_y   - self.bottomleft_y)/(self.topright_y-self.bottomleft_y);
        const CAMERA_HASH_SIZE_F: f32 = CAMERA_HASH_SIZE as f32;
        return (
            (clamp(left, 0.,0.999)*CAMERA_HASH_SIZE_F) as usize,
            (clamp(right,0.,0.999)*CAMERA_HASH_SIZE_F) as usize,
            (clamp(down ,0.,0.999)*CAMERA_HASH_SIZE_F) as usize,
            (clamp(up   ,0.,0.999)*CAMERA_HASH_SIZE_F) as usize
        );
    }
}
