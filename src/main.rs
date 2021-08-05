mod laplacian_branching;




fn main() {
    let mut sim = laplacian_branching::LaplacianBranchingSim::init(160, 100, 10.0, 5.0, 10.0);
    sim.image().save("test.png").ok();
    sim.grow();
    sim.image().save("test2.png").ok();
    sim.grow();
    sim.image().save("test3.png").ok();
    sim.grow();
    sim.image().save("test4.png").ok();
    sim.grow();
    sim.image().save("test5.png").ok();
    sim.grow();
    sim.image().save("test6.png").ok();
}
