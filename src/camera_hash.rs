
use crate::utils::INF;
use crate::bounding_box::BoundingBox;

const MAX_OBJECTS_PER_CELL: usize = 16;

#[derive(Copy,Clone)]
pub struct CellList{
    arr: [usize;MAX_OBJECTS_PER_CELL],
    count: usize,
}

impl CellList {
    pub fn new() -> Self { Self{arr: [0;MAX_OBJECTS_PER_CELL],count: 0} }
    pub fn add(&mut self,idx: usize) -> () {
        self.arr[self.count] = idx;
        self.count+=1;
    }
}

const CAMERA_HASH_SIZE: usize = 100;
pub struct CameraHash{
    cells: [[CellList;CAMERA_HASH_SIZE];CAMERA_HASH_SIZE],
    bottomleft_x: f32,
    bottomleft_y: f32,
    topright_x: f32,
    topright_y: f32,
}

impl CameraHash{
    pub fn new() -> Self {
        Self{
            cells: [[CellList::new();CAMERA_HASH_SIZE];CAMERA_HASH_SIZE],
            bottomleft_x: INF,
            bottomleft_y: INF,
            topright_x: -INF,
            topright_y: -INF,
        }
    }
    pub fn set_size(&mut self,bb: &BoundingBox) -> (){
        self.bottomleft_x = bb.bottomleft_x.min(self.bottomleft_x);
        self.bottomleft_y = bb.bottomleft_y.min(self.bottomleft_y);
        self.topright_x   = bb.topright_x.max(self.topright_x);
        self.topright_y   = bb.topright_y.max(self.topright_y);
    }
    pub fn add(&mut self,obj_idx: usize,bb: &BoundingBox) -> (){
        let x_step = (self.topright_x-self.bottomleft_x)/(CAMERA_HASH_SIZE as f32);
        let y_step = (self.topright_y-self.bottomleft_y)/(CAMERA_HASH_SIZE as f32);
        let min_i = (bb.bottomleft_x/x_step) as usize;
        let max_i = ((bb.topright_x/x_step) as usize).min(CAMERA_HASH_SIZE-1);
        let min_j = (bb.bottomleft_y/y_step) as usize;
        let max_j = ((bb.topright_y/y_step) as usize).min(CAMERA_HASH_SIZE-1);
        for i in min_i..=max_i{
            for j in min_j..=max_j{
                self.cells[i][j].add(obj_idx);
            }
        }
    }
}