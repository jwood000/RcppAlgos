#pragma once

#include <cstdint>
#include <vector>
#include <gmpxx.h>
#include <array>

// ************************** Count Function Examples *************************
//
// Here are the corresponding count functions that one would use for each
// example below:
//
// NotPartition   :
// RepStdAll      :  CountPartRep(20)
// RepNoZero      :  CountPartRepLen(20, 5)
// RepShort       :  CountPartRepLen(23, 3)
// RepCapped      :  CountPartRepLenCap(14, 3, 10) N.B. Get first part:
//                    (3, 5, 12); Map to match(c(3, 5, 12), 3:12) -->>
//                    (1, 3, 10) -->> sum(c(1, 3, 10)) = 14
//
// DstctStdAll    :  CountPartDistinct(20)
// DstctMultiZero :  CountPartsDistinctMultiZero(c(0, 0, 1, 2, 17), 20, 5)
// DstctOneZero   :  CountPartDistinctLen(25, 5) N.B. Add 1 to each element to
//                     obtain new target = 25
//
// DstctNoZero    :  CountPartDistinctLen(20, 5)
// DstctCapped    :  CountPartDistinctLenCap(20, 4, 9)
// DstctCappedMZ  :  CountPartsDistinctMultiZeroCap(c(0, 0, 9, 11), 20, 4, 11)
// LengthOne      :  1 or 0
// Multiset       :  CountPartsMultiset(rep(1:3, 5), c(1, 2, 2, 15))
// CompMultiset   :  CountPartsMultiset(rep(1:3, 5), c(1, 2, 2, 15), true)
//
// ****************************************************************************
//
//
// ************************** Definitions w/ Examples *************************
//
// RepStdAll      : Get all partitions. E.g. tar = 20;
//                   startZ = c(0, 0, 0, 0, 20): CountPartRep(20)
//
// RepNoZero      : E.g. tar = 20; startZ = c(1, 1, 1, 1, 15):
//                   CountPartRepLen(20, 5)
//
// RepShort       : Case where width isn't maximized E.g. tar = 20;
//                   startZ = c(0, 0, 20)
//
// RepCapped      : E.g. tar = 20 of width = 3 from the integers 3:12:
//                   CountPartRepCap(14, 3, 10)
//
// DstctStdAll    : Get all distinct partitions (0 can repeat) E.g. tar = 20;
//                   startZ = c(0, 0, 0, 0, 20)
//
// DstctMultiZero : Case where startZ doesn't maximize 0's. E.g.
//                   tar = 20 startZ = c(0, 0, 1, 2, 17)
//
// DstctOneZero   : Similar to above but can occur when IsMult = FALSE.
//                   E.g. tar = 20 startZ = c(0, 1, 2, 3, 14)
//
// DstctNoZero    : E.g. tar = 20 startZ = c(1, 2, 3, 4, 10)
//
// DstctCapped    : E.g. tar = 20, m = 4, from 1:9 gives startZ = c(1, 2, 8, 9)
//
// DstctCappedMZ. : E.g. tar = 20, m = 4, from 0:11, freqs = c(2, rep(1, 11))
//                   gives startZ = c(0, 0, 9, 11)
//
// LengthOne      : Any partition/composition when m = 1
//
// Multiset       : Partitions of non-trivial multisets. Non-trivial here means
//                   elements other than 0 have multiplicity > 1
//
// CoarseGrained  : This is equal to ConstraintType::PartitionEsque
//
// CompRepNoZero  : These are standard compositions with repetition.
//
// CompRepWeak    : Same as above but we allow terms of the seq to be zero.
//
// CmpRpZroNotWk  : This one is a little tricky. We have compositions, however
//                   we only want to see permutations of non-zero values.
//
// CmpDstctNoZero : Standard compositions with distinct parts:
//                   E.g. tar = 20; m = 5; startZ = c(1, 2, 3, 4, 10)
//
// CmpDstctZNotWk : Standard compositions with distinct parts and one or more
//                   zeros. Only non-zero values are considered when
//                   determining the next iteration.
//
// CmpDstctMZWeak : Same as above however we allow terms to be zero.
//
// CompMultiset   : Compositions of non-trivial multisets. Non-trivial here
//                   means elements other than 0 have multiplicity > 1
//
// NotPartition
//
// ****************************************************************************

enum class PartitionType {
    RepStdAll      = 0,
    RepNoZero      = 1,
    RepShort       = 2,
    RepCapped      = 3,
    DstctStdAll    = 4,
    DstctMultiZero = 5,
    DstctOneZero   = 6,
    DstctNoZero    = 7,
    DstctCapped    = 8,
    DstctCappedMZ  = 9,
    LengthOne      = 10,
    Multiset       = 11,
    CoarseGrained  = 12,
    CompRepNoZero  = 13,
    CompRepWeak    = 14,
    CmpRpZroNotWk  = 15,
    CmpDstctNoZero = 16,
    CmpDstctZNotWk = 17,
    CmpDstctMZWeak = 18,
    CompMultiset   = 19,
    NotPartition   = 20
};

constexpr const char* PrintPTypes[] = {
    "RepStdAll",
    "RepNoZero",
    "RepShort",
    "RepCapped",
    "DstctStdAll",
    "DstctMultiZero",
    "DstctOneZero",
    "DstctNoZero",
    "DstctCapped",
    "DstctCappedMZ",
    "LengthOne",
    "Multiset",
    "CoarseGrained",
    "CompRepNoZero",
    "CompRepWeak",
    "CmpRpZroNotWk",
    "CmpDstctNoZero",
    "CmpDstctZNotWk",
    "CmpDstctMZWeak",
    "CompMultiset",
    "NotPartition"
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
    mpz_class bigCount;
    bool isGmp = false;
    bool isRep = false;
    bool isMult = false;
    bool isDist = false;
    bool isComb = false;
    bool isPart = false;
    bool isComp = false;      // Are we dealing with compositions?
    bool isWeak = false;      // Do we allow terms of the sequence to be zero?
                              //
                              //     See: https://en.wikipedia.org/wiki/Composition_(combinatorics)
                              //
    bool allOne = false;      // When we have multisets with the pattern:
                              //
                              //     freqs = c(n, rep(1, p))
                              //
                              // This reduces to distinct
                              // partitions/compositions of differing widths.
                              //
                              // allOne translates to:
                              // "Every multiplicity is one expect the first element"
                              //
    bool mIsNull = false;     // Is the width provided by the user
    bool solnExist = false;   //
    bool includeZero = false; // Is the leading element zero?
    bool mapIncZero = false;  // There are some cases where includeZero = true,
                              // however after mapping, we don't have any
                              // zeros. E.g. tar = 20; startZ = c(0, 0, 0, 20);
                              // repetition = TRUE -->> mapTar = 24
    bool numUnknown = true;
    std::vector<int> startZ;
    std::int64_t cap = 0;
    std::int64_t shift = 0;
    std::int64_t slope = 0;
    std::int64_t target = 0;
    PartitionType ptype = PartitionType::NotPartition;
};
