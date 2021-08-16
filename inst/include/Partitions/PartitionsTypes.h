#ifndef PARTITION_TYPES_H
#define PARTITION_TYPES_H

#include <gmp.h>
#include <array>
#include <vector>
#include <cstdint>

// Here are the corresponding functions that one would use for each example below:
//
// NotPartition   :
// RepStdAll      :  CountPartRep(20)
// RepNoZero      :  CountPartRepLen(20, 5)
// RepShort       :  CountPartRepLen(23, 3)
// RepCapped      :  CountPartRepLenCap(14, 3, 10) N.B. Get first part: (3, 5, 12); Map to match(c(3, 5, 12), 3:12)
// DstctStdAll    :  CountPartDistinct(20)
// DstctMultiZero :  CountPartsDistinctMultiZero(c(0, 0, 1, 2, 17), 20, 5)
// DstctOneZero   :  CountPartDistinctLen(25, 5) N.B. Add 1 to each element to obtain new target = 25
// DstctNoZero    :  CountPartDistinctLen(20, 5)
// DstctCapped    :  CountPartDistinctLenCap(20, 4, 9)
// DstctCappedMZ  :  CountPartsDistinctMultiZeroCap(c(0, 0, 9, 11), 20, 4, 11)
// Multiset       :  CountPartMultiset(rep(1:3, 6), c(1, 2, 2, 15), 4, 20 - 1, 4 - 1)

enum class PartitionType {
    RepStdAll,       // Get all partitions. E.g. tar = 20 startZ = c(0, 0, 0, 0, 20): CountPartRep(20)
    RepNoZero,       // E.g. tar = 20 startZ = c(1, 1, 1, 1, 15): CountPartRepLen(20, 5)
    RepShort,        // Case where width isn't maximized E.g. tar = 20 startZ = c(0, 0, 20)
    RepCapped,       // E.g. tar = 20 of width = 3 from the integers 3:12: CountPartRepCap(14, 3, 10)
    DstctStdAll,     // Get all distinct partitions (0 can repeat) E.g. tar = 20 startZ = c(0, 0, 0, 0, 20)
    DstctMultiZero,  // Case where startZ doesn't maximize 0's. E.g. tar = 20 startZ = c(0, 0, 1, 2, 17)
    DstctOneZero,    // Similar to above but can occur when IsMult = FALSE. E.g. tar = 20 startZ = c(0, 1, 2, 3, 14)
    DstctNoZero,     // E.g. tar = 20 startZ = c(1, 2, 3, 4, 10)
    DstctCapped,     // E.g. tar = 20, m = 4, from 1:9 gives startZ = c(1, 2, 8, 9)
    DstctCappedMZ,   // E.g. tar = 20, m = 4, from 0:11, freqs = c(2, rep(1, 11)) gives startZ = c(0, 0, 9, 11)
    Multiset,        // Partitions of non-trivial multisets
    CoarseGrained,   // This is equivalent to ConstraintType::PartitionEsque
    NotPartition
};

const std::array<PartitionType, 3> RepPTypeArr{{
    PartitionType::RepStdAll, PartitionType::RepNoZero, PartitionType::RepShort
}};

const std::array<PartitionType, 5> DistPTypeArr{{
    PartitionType::DstctStdAll,  PartitionType::DstctMultiZero,
    PartitionType::DstctOneZero, PartitionType::DstctNoZero,
    PartitionType::DstctCappedMZ
}};

struct PartDesign {
    int width = 0;
    int mapTar = 0; // mapped target value
    double count = 0;
    mpz_t bigCount;
    bool isGmp = false;
    bool isRep = false;
    bool isMult = false;
    bool isPart = false;
    bool allOne = false;
    bool mIsNull = false;
    bool solnExist = false;
    bool includeZero = false;
    bool mapIncZero = false;
    bool numUnknown = false;
    std::vector<int> startZ;
    std::int64_t shift = 0;
    std::int64_t slope = 0;
    std::int64_t target = 0;
    PartitionType ptype = PartitionType::NotPartition;
};

#endif
