mod laplacian_branching;




fn main() {
    let mut sim = laplacian_branching::LaplacianBranchingSim::init(320, 240, 10.0, 5.0, 10.0);
    for i in 1..8 {
        sim.image().save(format!("sim_{}.png", i)).ok();
        sim.grow()
    }
}
