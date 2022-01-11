
use crate::utils::INF;
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
        let x_step = (self.topright_x-self.bottomleft_x)/(CAMERA_HASH_SIZE as f32);
        let y_step = (self.topright_y-self.bottomleft_y)/(CAMERA_HASH_SIZE as f32);
        let left_dis  = bb.bottomleft_x - self.bottomleft_x;
        let right_dis =   bb.topright_x - self.bottomleft_x;
        let down_dis  = bb.bottomleft_y - self.bottomleft_y;
        let up_dis    =   bb.topright_y - self.bottomleft_y;
        let min_i =  (left_dis/x_step).max(0.) as usize;
        let max_i = ((right_dis/x_step) as usize).min(CAMERA_HASH_SIZE-1);
        let min_j =  (down_dis/y_step).max(0.) as usize;
        let max_j = ((up_dis/y_step) as usize).min(CAMERA_HASH_SIZE-1);
        return (min_i,max_i,min_j,max_j);
    }
}
