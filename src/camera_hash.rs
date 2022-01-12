
use crate::utils::clamp;
use crate::bounding_box::BoundingBox;

pub trait CellEmptyInitializable: Copy {
    fn empty(&mut self) -> ();
}

const CAMERA_HASH_SIZE: usize = 100;
#[derive(Debug)]
pub struct CameraHash<Cell>{
    pub cells: [[Cell;CAMERA_HASH_SIZE];CAMERA_HASH_SIZE],
    pub bb: BoundingBox,
}

impl <Cell: CellEmptyInitializable> CameraHash<Cell>{
    pub fn empty(&mut self){
        for i in 0..CAMERA_HASH_SIZE{
            for j in 0..CAMERA_HASH_SIZE{
                self.cells[i][j].empty();
            }
        }
        self.bb = BoundingBox::draw_always();
    }
    pub fn get_indexes(&self,bb: &BoundingBox) -> (usize,usize,usize,usize){
        let width_inv  = 1./(self.bb.topright_x-self.bb.bottomleft_x);
        let height_inv = 1./(self.bb.topright_y-self.bb.bottomleft_y);
        let left  = (bb.bottomleft_x - self.bb.bottomleft_x)*width_inv;
        let right = (bb.topright_x   - self.bb.bottomleft_x)*width_inv;
        let down  = (bb.bottomleft_y - self.bb.bottomleft_y)*height_inv;
        let up    = (bb.topright_y   - self.bb.bottomleft_y)*height_inv;
        const CAMERA_HASH_SIZE_F: f32 = CAMERA_HASH_SIZE as f32;
        let left  = (clamp(left, 0.,0.999)*CAMERA_HASH_SIZE_F) as usize;
        let right = (clamp(right,0.,0.999)*CAMERA_HASH_SIZE_F) as usize;
        let down  = (clamp(down ,0.,0.999)*CAMERA_HASH_SIZE_F) as usize;
        let up    = (clamp(up   ,0.,0.999)*CAMERA_HASH_SIZE_F) as usize;
        return (left.min(right),left.max(right),down.min(up),down.max(up));
    }
}
