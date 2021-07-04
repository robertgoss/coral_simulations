// Simulate a branching coral by assuming nutirent 
// suspension follows a laplacian distribution with 
// the substraint and exisiting structure absorbing nutrients 
// supplied from the water layer.
//
// Each branch is assumed to grow towards the highest concentation
// of nutriants and split if the concentration is even.

use ndarray::{Array, Ix2, ShapeBuilder};

// Skeleton of coral
struct Node {
    point : (i32, i32),
    children : Vec<Node>
}
struct Skeleton {
    root : Node,
    thickness : u32
}

// Nutrient concentraion
struct Concentraion {
    phi : Array<f64, Ix2>,
    x_offset : u32
}


impl Concentraion {
    fn background(width: u32, height: u32) -> Concentraion {
        Concentraion {
            phi : Array::from_shape_fn(
                (width as usize, height as usize).f(),
                |(i,_)| f64::from(i as u32) / f64::from(height)
            ),
            x_offset : width / 2
        }
    }

    fn init(width: u32, height: u32, coral : &Skeleton) -> Concentraion {
        let mut concentraion = Concentraion::background(width, height);

        concentraion
    }
}


impl Skeleton {
    // Create initial skeleton with one piece 
    pub fn init(initial_height : i32, thickness : u32) -> Skeleton {
        Skeleton {
            root : Node {
                point : (0, 0),
                children : vec!(
                    Node { point : (0, initial_height), children : Vec::new()}
                )
            },
            thickness : thickness
        }
    }
}