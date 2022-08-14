#ifndef PARTITION_TYPES_H
#define PARTITION_TYPES_H

#include <cstdint>
#include <vector>
#include <gmp.h>
#include <array>

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
// LengthOne      :  1 or 0
// Multiset       :  CountPartMultiset(rep(1:3, 6), c(1, 2, 2, 15), 4, 20 - 1, 4 - 1)

enum class PartitionType {
    RepStdAll      = 0,  // Get all partitions. E.g. tar = 20 startZ = c(0, 0, 0, 0, 20): CountPartRep(20)
    RepNoZero      = 1,  // E.g. tar = 20 startZ = c(1, 1, 1, 1, 15): CountPartRepLen(20, 5)
    RepShort       = 2,  // Case where width isn't maximized E.g. tar = 20 startZ = c(0, 0, 20)
    RepCapped      = 3,  // E.g. tar = 20 of width = 3 from the integers 3:12: CountPartRepCap(14, 3, 10)
    DstctStdAll    = 4,  // Get all distinct partitions (0 can repeat) E.g. tar = 20 startZ = c(0, 0, 0, 0, 20)
    DstctMultiZero = 5,  // Case where startZ doesn't maximize 0's. E.g. tar = 20 startZ = c(0, 0, 1, 2, 17)
    DstctOneZero   = 6,  // Similar to above but can occur when IsMult = FALSE. E.g. tar = 20 startZ = c(0, 1, 2, 3, 14)
    DstctNoZero    = 7,  // E.g. tar = 20 startZ = c(1, 2, 3, 4, 10)
    DstctCapped    = 8,  // E.g. tar = 20, m = 4, from 1:9 gives startZ = c(1, 2, 8, 9)
    DstctCappedMZ  = 9,  // E.g. tar = 20, m = 4, from 0:11, freqs = c(2, rep(1, 11)) gives startZ = c(0, 0, 9, 11)
    LengthOne      = 10, // Any partition when m = 1
    Multiset       = 11, // Partitions of non-trivial multisets
    CoarseGrained  = 12, // This is equivalent to ConstraintType::PartitionEsque
    NotPartition   = 13
};

const std::array<PartitionType, 4> RepPTypeArr{{
    PartitionType::RepStdAll, PartitionType::RepNoZero,
    PartitionType::RepShort, PartitionType::RepCapped
}};

const std::array<PartitionType, 6> DistPTypeArr{{
    PartitionType::DstctStdAll,  PartitionType::DstctMultiZero,
    PartitionType::DstctOneZero, PartitionType::DstctNoZero,
    PartitionType::DstctCapped, PartitionType::DstctCappedMZ
}};

struct PartDesign {
    int width = 0;
    int mapTar = 0; // mapped target value
    double count = 0;
    mpz_t bigCount;
    bool isGmp = false;
    bool isRep = false;
    bool isMult = false;
    bool isComb = false;
    bool isPart = false;
    bool isComp = false;      // Are we dealing with compositions?
    bool isWeak = false;      // Do we allow terms of the sequence to be zero?
                              //  See: https://en.wikipedia.org/wiki/Composition_(combinatorics)
    bool allOne = false;      // When we have multisets with the pattern:
                              //     freqs = c(n, rep(1, p))
                              // This reduces to distinct
                              // partitions/compositions of differing widths.
                              // allOne translates to:
                              // "Every multiplicity is one expect the first element"

    bool mIsNull = false;     // Is the width provided by the user
    bool solnExist = false;
    bool includeZero = false; // Is the leading element zero?
    bool mapIncZero = false;
    bool numUnknown = true;
    std::vector<int> startZ;
    std::int64_t cap = 0;
    std::int64_t shift = 0;
    std::int64_t slope = 0;
    std::int64_t target = 0;
    PartitionType ptype = PartitionType::NotPartition;
};

#endif
