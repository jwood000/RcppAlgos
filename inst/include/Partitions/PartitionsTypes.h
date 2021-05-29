#ifndef PARTITION_TYPES_H
#define PARTITION_TYPES_H

#include <array>
#include <vector>
#include <cstdint>

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

// Here are the corresponding functions that one would use for each example below:
//
// NotPartition :
// RepStdAll    :  CountPartRep(20)
// RepNoZero    :  CountPartRepLen(20, 5)
// RepShort     :  CountPartRepLen(23, 3)
// RepCapped    :  CountPartRepLenCap(14, 3, 10) N.B. Get first part: (3, 5, 12); Map to match(c(3, 5, 12), 3:12)
// DstctStdAll  :  CountPartDistinct(20)
// DstctShort   :  GetSpecialCount(c(0, 0, 20), 20, 3) N.B. We can't use the "add 1 trick" as zero is repeated.
// DstctSpecial :  GetSpecialCount(c(0, 0, 1, 2, 17), 20, 5) ... This would give c(1, 1, 21) which isn't distinct
// DstctOneZero :  CountPartDistinctLen(25, 5) N.B. Add 1 to each element to obtain new target = 25
// DstctNoZero  :  CountPartDistinctLen(20, 5)
// DistCapped   :  CountPartDistinctLenCap(20, 4, 9)
// Multiset     :  CountPartMultiset(rep(1:3, 6), c(1, 2, 2, 15), 4, 20 - 1, 4 - 1)

enum class PartitionType {
    RepStdAll,    // Get all partitions. E.g. tar = 20 startZ = c(0, 0, 0, 0, 20): CountPartRep(20)
    RepNoZero,    // E.g. tar = 20 startZ = c(1, 1, 1, 1, 15): CountPartRepLen(20, 5)
    RepShort,     // Case where width isn't maximized E.g. tar = 20 startZ = c(0, 0, 20)
    RepCapped,    // E.g. tar = 20 of width = 3 from the integers 3:12: CountPartRepCap(14, 3, 10)
    DstctStdAll,  // Get all distinct partitions (0 can repeat) E.g. tar = 20 startZ = c(0, 0, 0, 0, 20)
    DstctShort,   // Case where width isn't maximized E.g. tar = 20 startZ = c(0, 0, 20)
    DstctSpecial, // Case where startZ doesn't maximize 0's. E.g. tar = 20 startZ = c(0, 0, 1, 2, 17)
    DstctOneZero, // Similar to above but can occur when IsMult = FALSE. E.g. tar = 20 startZ = c(0, 1, 2, 3, 14)
    DstctNoZero,  // E.g. tar = 20 startZ = c(1, 2, 3, 4, 10)
    DistCapped,   // E.g. tar = 20, m = 4, from 1:9 gives startZ = c(1, 2, 8, 9)
    Multiset,     // Partitions of non-trivial multisets
    NotPartition
};

const std::array<PartitionType, 3> RepPTypeArr{{
    PartitionType::RepStdAll, PartitionType::RepNoZero, PartitionType::RepShort
}};

const std::array<PartitionType, 5> DistPTypeArr{{
    PartitionType::DstctStdAll,  PartitionType::DstctShort,
    PartitionType::DstctSpecial, PartitionType::DstctOneZero,
    PartitionType::DstctNoZero
}};

struct PartDesign {
    int width = 0;
    int mapTar = 0; // mapped target value
    double count = 0;
    bool isRep = false;
    bool isMult = false;
    bool allOne = false;
    bool mIsNull = false;
    bool solnExist = false;
    bool mapZeroFirst = false;
    std::vector<int> startZ;
    std::int64_t shift = 0;
    std::int64_t slope = 0;
    std::int64_t target = 0;
    ConstraintType ctype = ConstraintType::General;
    PartitionType ptype = PartitionType::NotPartition;
};

#endif
