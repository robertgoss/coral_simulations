mod laplacian_branching;




fn main() {
    let skeleton = laplacian_branching::Skeleton::init(5.0, 2.5);
    let mut conc = skeleton.get_concentration(80, 50);
    conc.image().save("test.png").unwrap();
    conc.diffuse();
}
