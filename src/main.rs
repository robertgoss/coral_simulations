mod laplacian_branching;




fn main() {
    let skeleton = laplacian_branching::Skeleton::init(5.0, 2.5);
    let conc = skeleton.get_concentration(800, 500);
    conc.image().save("test.png").unwrap();
}
