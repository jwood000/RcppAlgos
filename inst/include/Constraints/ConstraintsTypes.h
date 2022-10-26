#pragma once

// SpecialCnstrnt = 1: For this case, we have either lower = true,
// which means the user is looking for combinations/permutations in a
// particular range or if we have constraintFun = "prod" and negative values
// are present. When the latter occurs, it is non-trivial to produce a loose
// monotonic sequence of results.

// ConstraintType - General: Cannot be reduced to a specific partiion case.
// The most general algorithm is used in this case.

// PartitionEsque = 3: Can't be reduced to an integer partition but still
// has similarities to the more general subset sum problem. E.g.
// v = rnorm(20, mean = 10.5), m = 4, rep = TRUE, tar = c(11.005, 11.15),
//  comparisonFun = c(">", "<"), constraintFun = "mean"

// PartMapping = 4: Occurs when non-standard input can be reduced to a general
// integer partition: E.g. v = seq(200, 300, 5), tar = 1100, m = 4, rep = TRUE

// PartStandard = 5: The algorithms used to generate the base structures for
// PartMapping and PartStandard are the same, however PartMapping requires
// additional mappings. E.g. For the example given above, the analog for this
// case would be: v = 0:20, m = 4, tar = 60, rep = TRUE.
//
// The 3 examples below are isomorphically equivalent:
//
// ## PartStandard = 5
// comboGeneral(0:20, 4, T, constraintFun = "sum",
//                          comparisonFun = "==",
//                          limitConstraints = 60)
//
// ## This is the same as above only there is no zero
// comboGeneral(21, 4, T, constraintFun = "sum",
//                          comparisonFun = "==",
//                          limitConstraints = 64)
//
// ## PartMapping = 4
// comboGeneral(seq(200, 300, 5), 4, T, constraintFun = "sum",
//                                      comparisonFun = "==",
//                                      limitConstraints = 1100)

enum class ConstraintType {
    NoConstraint   = 0,
    SpecialCnstrnt = 1,
    General        = 2,
    PartitionEsque = 3,
    PartMapping    = 4,
    PartStandard   = 5
};

enum class FunType {
    Min  = 1,
    Max  = 2,
    Sum  = 3,
    Prod = 4,
    Mean = 5
};
