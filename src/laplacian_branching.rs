// Simulate a branching coral by assuming nutirent 
// suspension follows a laplacian distribution with 
// the substraint and exisiting structure absorbing nutrients 
// supplied from the water layer.
//
// Each branch is assumed to grow towards the highest concentation
// of nutriants and split if the concentration is even.

use ndarray::{Array, Ix2, ShapeBuilder,s , Zip};
use image::{RgbImage, Rgb};

use geo::{Coordinate, Line};
use geo::algorithm::euclidean_distance::EuclideanDistance;

type Point2D = Coordinate<f64>;

// Skeleton of coral
struct Node {
    point : Point2D,
    children : Vec<Node>
}
pub struct Skeleton {
    root : Node,
    thickness : f64
}

// Nutrient concentraion
pub struct Concentraion {
    phi : Array<f64, Ix2>,
    x_offset : u32
}


impl Concentraion {

    fn init(width: u32, height: u32, coral : &Skeleton) -> Concentraion {
        let x_offset = width / 2;
        Concentraion {
            phi : Array::from_shape_fn(
                (width as usize, height as usize).f(),
                |(i,j)| if coral.contains((i as f64) - (x_offset as f64) , (height as f64) - (j as f64)) {
                            1.0
                        } else {
                            1.0 - f64::from(j as u32) / f64::from(height)
                        }
            ),
            x_offset : x_offset
        }
    }

    pub fn diffuse(self : &mut Self) {
        // Note can do better for now do fixed with small timestep
        for i in 1.. {
            self.image().save(format!("img/iter_{}.png", i)).ok();
            if self.step_laplacian(0.1) {
                break;
            }
        }
    }

    fn step_laplacian(self : &mut Self, step : f64) -> bool {
        let old_phi : Array<f64, Ix2> = self.phi.clone();
        let left = old_phi.slice(s![2.., 1..-1]);
        let right = old_phi.slice(s![0..-2, 1..-1]);
        let up = old_phi.slice(s![1..-1, 2..]);
        let down = old_phi.slice(s![1..-1, 0..-2]);
        let mut mid = self.phi.slice_mut(s![1..-1, 1..-1]);
        Zip::from(&mut mid).and(left).and(right).and(up).and(down).for_each(
            |m, &l, &r, &u, &d| {
                *m = if *m == 0.0 {
                    0.0 // Solid case
                } else {
                    *m + step * (l+r+u+d - 4.0 * (*m))
                }
            }
        );
        let diff = Zip::from(&old_phi).and(&self.phi).fold(0.0, |acc, a, b| acc + (a - b).abs());
        diff < step
    }

    pub fn image(self : &Self) -> RgbImage {
        let (width, height) = self.phi.dim();
        let mut img = RgbImage::new(width as u32, height as u32);
        for x in 0..width {
            for y in 0..height {
                let val : u8 = (255.0 * self.phi.get( (x,y) ).unwrap()) as u8;
                img.put_pixel(x as u32, y as u32, Rgb([val, 0, 0]));
            }
        }
        img
    }
}


impl Skeleton {
    // Create initial skeleton with one piece 
    pub fn init(initial_height : f64, thickness : f64) -> Skeleton {
        Skeleton {
            root : Node {
                point : Point2D { x: 0.0, y: 0.0},
                children : vec!(
                    Node { point : Point2D {x :0.0, y : initial_height}, children : Vec::new()}
                )
            },
            thickness : thickness
        }
    }

    pub fn get_concentration(self : &Self, width: u32, height: u32) -> Concentraion {
        Concentraion::init(width, height, self)
    }

    fn contains(self : &Self, x : f64, y : f64) -> bool {
        self.root.point_within(&Point2D{x: x, y: y}, self.thickness)
    }
}

impl Node {
    fn point_within(self : &Self, point : &Point2D, dist : f64) -> bool {
        self.point_within_direct(point, dist) ||
          self.children.iter().any(
              |child| child.point_within(point, dist)
          )
    }

    fn point_within_direct(self : &Self, point : &Point2D, dist : f64) -> bool {
        self.children.iter().map(
            |child| Line::new(self.point, child.point)
        ).any(
            |branch| point.euclidean_distance(&branch) < dist
        )
    }
}