// Simulate a branching coral by assuming nutirent 
// suspension follows a laplacian distribution with 
// the substraint and exisiting structure absorbing nutrients 
// supplied from the water layer.
//
// Each branch is assumed to grow towards the highest concentation
// of nutriants and split if the concentration is even.

use ndarray::{Array, Ix2, ShapeBuilder,s , Zip};
use image::{RgbImage, Rgb};
use imageproc::pixelops;
use imageproc::drawing::{draw_antialiased_line_segment_mut, draw_filled_circle_mut};

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
    thickness : f64,
    growth : f64
}

// Nutrient concentraion
pub struct Concentration {
    phi : Array<f64, Ix2>,
    rad : f64,
    x_offset : usize
}

// Simulation parameters
pub struct LaplacianBranchingSim {
    width: u32, 
    height: u32, 
    grad_rad : f64,
    skeleton : Skeleton,
    concentration : Concentration
}


impl Concentration {

    fn init(width: u32, height: u32, grad_rad : f64, coral : &Skeleton) -> Concentration {
        let x_offset = width / 2;
        Concentration {
            phi : Array::from_shape_fn(
                (width as usize, height as usize).f(),
                |(i,j)| if coral.contains((i as f64) - (x_offset as f64) , j as f64) {
                            0.0
                        } else {
                            f64::from(j as u32) / f64::from(height)
                        }
            ),
            rad : grad_rad,
            x_offset : x_offset as usize
        }
    }

    pub fn diffuse(self : &mut Self) {
        // Note can do better for now do fixed with small timestep
        for _ in 1.. {
            if self.step_laplacian(0.1) {
                break;
            }
        }
    }

    fn step_laplacian(self : &mut Self, step : f64) -> bool {
        let old_phi : Array<f64, Ix2> = self.phi.clone();
        // Do discrete time step in the middle
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
        // Special handling of sides for periodic sides
        let left_side = old_phi.slice(s![-2, 1..-1]);
        let right_side = old_phi.slice(s![1, 1..-1]);
        let up_side= old_phi.slice(s![0, 2..]);
        let down_side = old_phi.slice(s![0, 0..-2]);
        let mut side1 = self.phi.slice_mut(s![0, 1..-1]);
        Zip::from(&mut side1).and(left_side).and(right_side).and(up_side).and(down_side).for_each(
            |m, &l, &r, &u, &d| { *m = *m + step * (l+r+u+d - 4.0 * (*m)) }
        );
        let mut side2 = self.phi.slice_mut(s![0, 1..-1]);
        Zip::from(&mut side2).and(left_side).and(right_side).and(up_side).and(down_side).for_each(
            |m, &l, &r, &u, &d| { *m = *m + step * (l+r+u+d - 4.0 * (*m)) }
        );
        let diff = Zip::from(&old_phi).and(&self.phi).fold(0.0, |acc, a, b| acc + (a - b).abs());
        diff < step
    }

    fn concentration(self : &Self, x : i64, y : i64) -> f64 {
        *self.phi.get( (x as usize, y as usize) ).unwrap_or(&0.0)
    }

    fn average_gradient(self : &Self, point : &Point2D) -> (f64, f64) {
        // Get box range
        let x_min = (point.x + self.x_offset as f64 - self.rad).floor() as i64;
        let x_max = (point.x + self.x_offset as f64 + self.rad).ceil() as i64;
        let y_min = (point.y - self.rad).floor() as i64;
        let y_max = (point.y + self.rad).ceil() as i64;
        let local_point = Point2D {
            x: point.x + self.x_offset as f64,
            y: point.y
        };
        // Accumalate the gradients - this was easier as a loop - looking for a better way.
        let mut count = 0;
        let mut dx = 0.0;
        let mut dy = 0.0;
        for x in x_min..=x_max {
            for y in y_min..=y_max {
                if self.concentration(x,y) != 0.0 {
                    if local_point.euclidean_distance(&Point2D { x: x as f64, y: y as f64 }) < self.rad {
                        dx += self.concentration(x+1,y) - self.concentration(x-1,y);
                        dy += self.concentration(x,y+1) - self.concentration(x,y-1);
                        count += 1;
                    }
                }
            }
        }
        (dx / (count as f64), dy / (count as f64))
    }

    fn interpolated_concentration(self : &Self, point : &Point2D) -> f64 {
        let x = (point.x + self.x_offset as f64).floor() as i64;
        let y = (point.y).floor() as i64;
        let c11 = self.concentration(x,y);
        let c12 = self.concentration(x,y+1);
        let c21 = self.concentration(x+1,y);
        let c22 = self.concentration(x+1,y+1);
        let u = point.x - point.x.floor();
        let v = point.y - point.y.floor();
        let c1 = (1.0-v)*c11 + v*c12;
        let c2 = (1.0-v)*c21 + v*c22;
        (1.0-u)*c1 + u*c2
    }

    pub fn image(self : &Self) -> RgbImage {
        let (width, height) = self.phi.dim();
        let mut img = RgbImage::new(width as u32, height as u32);
        for x in 0..width {
            for y in 0..height {
                let val : u8 = (255.0 * self.phi.get( (x,y) ).unwrap()) as u8;
                img.put_pixel(x as u32, (height-1 - y) as u32, Rgb([val, 0, 0]));
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
            growth : initial_height,
            thickness : thickness
        }
    }

    fn grow(self : &mut Self, concentration : &Concentration) {
        self.root.grow(concentration, self.growth)
    }

    fn contains(self : &Self, x : f64, y : f64) -> bool {
        self.root.point_within(&Point2D{x: x, y: y}, self.thickness)
    }

    fn draw(self : &Self, image : &mut RgbImage) {
        self.root.draw(image);
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

    fn draw(self : &Self, image : &mut RgbImage) {
        let pixel = self.pixel(image.width(), image.height());
        for child in &self.children {
            draw_antialiased_line_segment_mut(
                image, 
                pixel, 
                child.pixel(image.width(), image.height()),
                Rgb([255, 255, 255]),
                pixelops::interpolate
            );
            child.draw(image);
        }
        draw_filled_circle_mut(image, pixel, 2, Rgb([0, 128, 128]));
    }

    fn pixel(self : &Self, width : u32, height : u32) -> (i32, i32) {
        let x = (self.point.x.round() as i32) + (width as i32 / 2);
        let y = (height as i32) - (self.point.y.round() as i32);
        (x, y)
    }

    fn grow_direct(self : &mut Self, length : f64, direction : &(f64,f64)) {
        let mag = ((direction.0 * direction.0) + (direction.1 * direction.1)).sqrt();
        if mag > 0.0 {
            let point = Point2D {
                x : self.point.x + (length / mag) * direction.0, 
                y : self.point.y + (length / mag) * direction.1
            };
            self.children.push(
                Node { point : point, children : Vec::new()}
            );
        }
    }

    fn grow(self : &mut Self, concentation : &Concentration, length : f64) {
        if self.children.is_empty() {
            let grad = concentation.average_gradient(&self.point);
            self.grow_direct(length, &grad);
        } else {
            self.children.iter_mut().for_each(
                |child| child.grow(concentation, length)
            );
        }
    } 
}

impl LaplacianBranchingSim {
    pub fn init(width: u32, height: u32, initial_height : f64, thickness : f64, grad_rad : f64) -> LaplacianBranchingSim {
        let skeleton = Skeleton::init(initial_height, thickness);
        let concentration = Concentration::init(width, height, grad_rad, &skeleton);
        LaplacianBranchingSim {
            width : width,
            height : height,
            grad_rad : grad_rad,
            skeleton : skeleton,
            concentration : concentration
        }
    }

    pub fn grow(self : &mut Self) {
        self.concentration.diffuse();
        self.skeleton.grow(&self.concentration);
        self.concentration = Concentration::init(
            self.width, 
            self.height, 
            self.grad_rad, 
            &self.skeleton
        );
    }

    pub fn image(self : &Self) -> RgbImage {
        let mut image : RgbImage = self.concentration.image();
        self.skeleton.draw(&mut image);
        image
    }
}