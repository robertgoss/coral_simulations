mod laplacian_branching;




fn main() {
    let sim = laplacian_branching::LaplacianBranchingSim::init(160, 100, 10.0, 5.0);
    sim.image().save("test.png").ok();
}
