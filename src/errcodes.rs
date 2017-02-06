#[derive(Copy,Clone,Debug)]
pub enum GeomErrorCode {
    TangentVectorsTooSmall,
    DegenerateOrSmallConic,
    VectorsNotInPlaneOfPoints,
    DegenerateCircle,
    JoinContinuityTooTight,
    BadKnots,
    MismatchedArraySizes,
    NotEnoughPointsForInterp,
    DuplicatePointData,
    InvalidPeriodicData,
    NotEnoughPointsForPeriodicInterp,
}
