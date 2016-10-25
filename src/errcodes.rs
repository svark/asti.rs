#[derive(Copy,Clone,Debug)]
pub enum GeomErrorCode {
    TangentVectorsTooSmall,
    DegenerateOrSmallConic,
    VectorsNotInPlaneOfPoints,
    DegenerateCircle,
    JoinContinuityTooTight,
}
