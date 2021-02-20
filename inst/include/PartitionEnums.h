#ifndef PARTITION_ENUMS_H
#define PARTITION_ENUMS_H

// PartitionEsque = 2: Can't be reduced to an integer partition but still has similarities 
// to the more general subset sum problem. E.g. v = rnorm(20, mean = 10.5), m = 4,
// rep = TRUE, tar = c(11.005, 11.15), comparisonFun = c(">", "<"), constraintFun = "mean"

// PartMapping = 3: Occurs when non-standard input can be reduced to a general integer
// partition: E.g. v = seq(200, 300, 5), tar = 1100, m = 4, rep = TRUE

// PartStandard = 4: The algorithms used to generate the base structures for PartMapping
// and PartStandard are the same, however PartMapping requires additional mappings. E.g.
// For the example given above, the analog for this case would be:
// v = 0:20, m = 4, tar = 60, rep = TRUE.
// 
// The 3 examples below are isomorphically equivalent:
//
// ## PartStandard = 4
// comboGeneral(0:20, 4, T, constraintFun = "sum",
//                          comparisonFun = "==",
//                          limitConstraints = 60)
//
// ## This is the same as above only there is no zero
// comboGeneral(21, 4, T, constraintFun = "sum",
//                          comparisonFun = "==",
//                          limitConstraints = 64)
//
// ## PartMapping = 3
// comboGeneral(seq(200, 300, 5), 4, T, constraintFun = "sum",
//                                      comparisonFun = "==",
//                                      limitConstraints = 1100)

enum class ConstraintType {
    General        = 1, // Cannot be reduced to a specific partiion case
    PartitionEsque = 2,
    PartMapping    = 3,
    PartStandard   = 4,
};

enum class PartitionType {
    Traditional   = 1, // Get all partitions. E.g. tar = 20 startZ = c(0, 0, 0, 0, 20)
    TradNoZero    = 2, // E.g. tar = 20 startZ = c(1, 1, 1, 1, 15)
    TradCapped    = 3, // E.g. tar = 20 from the integers 3:14
    DstctStdAll   = 4, // Get all distinct partitions (0 can repeat) E.g. tar = 20 startZ = c(0, 0, 0, 0, 20)
    DstctShort    = 5, // Case where startZ doesn't maximize width. E.g. tar = 20 startZ = c(0, 0, 20)
    DstctSpecial  = 6, // Case where startZ doesn't maximize 0's. E.g. tar = 20 startZ = c(0, 0, 1, 2, 17)
    DstctOneZero  = 7, // Similar to above but can occur when IsMult = FALSE. E.g. tar = 20 startZ = c(0, 1, 2, 3, 14)
    DstctNoZero   = 8, // E.g. tar = 20 startZ = c(1, 2, 3, 4, 10)
    DistCapped    = 9, // E.g. tar = 20, m = 4, from 1:9 gives startZ = c(1, 2, 8, 9)
};

enum class Sign {
    Positive = 1,
    Negitive = 2,
    MixedBag = 3,
};

#endif
