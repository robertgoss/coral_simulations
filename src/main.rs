mod laplacian_branching;




fn main() {
    let mut sim = laplacian_branching::LaplacianBranchingSim::init(480, 320, 10.0, 5.0, 10.0, 0.5);
    sim.image().save(format!("sim_{}.png", 0)).ok();
    for i in 1..32 {
        sim.grow();
        println!("Stage {} - {} growing tips", i, sim.growing_tips());
        sim.image().save(format!("sim_{}.png", i)).ok();

    }
}
