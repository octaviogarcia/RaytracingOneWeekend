
use crate::utils::clamp;
use crate::bounding_box::BoundingBox;

pub trait CellEmptyInitializable: Copy {
    fn empty(&mut self) -> ();
}

const CAMERA_HASH_SIZE: usize = 100;
const CAMERA_HASH_SIZE_F: f32 = CAMERA_HASH_SIZE as f32;
#[derive(Debug)]
pub struct CameraHash<Cell>{
    pub cells: [[Cell;CAMERA_HASH_SIZE];CAMERA_HASH_SIZE],
}

impl <Cell: CellEmptyInitializable> CameraHash<Cell>{
    pub fn empty(&mut self){
        for i in 0..CAMERA_HASH_SIZE{
            for j in 0..CAMERA_HASH_SIZE{
                self.cells[i][j].empty();
            }
        }
    }
    pub fn get_indexes(&self,bb: &BoundingBox) -> (usize,usize,usize,usize){//bb should be in range [0;1] from screen cords
        let left  = (clamp(bb.bottomleft_x, 0.,0.999)*CAMERA_HASH_SIZE_F) as usize;//@SPEED: Should I clamp?
        let right = (clamp(bb.topright_x,0.,0.999)*CAMERA_HASH_SIZE_F) as usize;
        let down  = (clamp(bb.bottomleft_y ,0.,0.999)*CAMERA_HASH_SIZE_F) as usize;
        let up    = (clamp(bb.topright_y,0.,0.999)*CAMERA_HASH_SIZE_F) as usize;
        return (left.min(right),left.max(right),down.min(up),down.max(up));//shouldn't really need to min() and max()
    }
    pub fn get_index(&self,u: f32,v: f32) -> (usize,usize){
        const CAMERA_HASH_SIZE_F: f32 = CAMERA_HASH_SIZE as f32;
        let i = (clamp(u,0.,0.999)*CAMERA_HASH_SIZE_F) as usize;
        let j = (clamp(v,0.,0.999)*CAMERA_HASH_SIZE_F) as usize;
        return (i,j);
    }
    pub fn at(&self,u: f32,v: f32) -> &Cell {
        let (i,j) = self.get_index(u,v);
        return &self.cells[i][j];
    }
}
