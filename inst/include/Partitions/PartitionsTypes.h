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
// RepStdAll      : CountPartRep(20)
// RepNoZero      : CountPartRepLen(20, 5)
// RepShort       : CountPartRepLen(23, 3)
// RepCapped      : CountPartRepLenCap(14, 3, 10) N.B. Get first part:
//                   (3, 5, 12); Map to match(c(3, 5, 12), 3:12) -->>
//                   (1, 3, 10) -->> sum(c(1, 3, 10)) = 14
//
// DstctStdAll    : CountPartDistinct(20)
// DstctMultiZero : CountPartsDistinctMultiZero(c(0, 0, 1, 2, 17), 20, 5)
// DstctOneZero   : CountPartDistinctLen(25, 5) N.B. Add 1 to each element to
//                   obtain new target = 25
//
// DstctNoZero    : CountPartDistinctLen(20, 5)
// DstctCapped    : CountPartDistinctLenCap(20, 4, 9)
// DstctCappedMZ  : CountPartsDistinctMultiZeroCap(c(0, 0, 9, 11), 20, 4, 11)
// LengthOne      : 1 or 0
// Multiset       : CountPartsMultiset(rep(1:3, 5), c(1, 2, 2, 15))
//
// CoarseGrained  : Currently, we utilize std::vector and push_back until we
//                   reach a terminating situation.
//
// CompRepNoZero  : compositionsCount(20, 5, TRUE) -->> CountCompsRepLen(20, 5)
// CompRepWeak    : compositionsCount(0:20, 5, TRUE, weak = TRUE) -->>
//                   CountCompsRepLen(25, 5) We add the length to the target
//
// CmpRpZroNotWk  : compositionsCount(0:20, 5, TRUE) -- >>
//                   CountCompsRepZNotWk(20, 5)
//
// CmpDstctNoZero : compositionsCount(20, 5) -->> CountCompsDistinctLen(20, 5)
// CmpDstctZNotWk : compositionsCount(0:20, 5, freqs = c(3, rep(1, 20))) -->>
//                   CountCompsDistinctMultiZero(20, 5, 0, 2)
//
// CmpDstctMZWeak : compositionsCount(0:20, 5, freqs = c(3, rep(1, 20)),
//                                    weak = TRUE) -->>
//                   CountCompsDistinctMZWeak(20, 5, 0, 2)
//
// CompMultiset   : CountPartsMultiset(rep(1:3, 5), c(1, 2, 2, 15), true)
//
// PrmRepPartNoZ  : permuteCount(
//                      1:20, 5, TRUE,
//                      constraintFun = "sum",
//                      comparisonFun = "==",
//                      limitConstraints = 20
//                  ) -->> CountCompsRepLen(20, 5)
//
//                  ## When zero is involved, we call the same compiled
//                  ## function just translated by the width
// PrmRepPart     : permuteCount(
//                      0:20, 5, TRUE,
//                      constraintFun = "sum",
//                      comparisonFun = "==",
//                      limitConstraints = 20
//                  ) -->> CountCompsRepLen(25, 5)
//
// PrmRepCapped   : permuteCount(
//                      1:7, 5, TRUE,
//                      constraintFun = "sum",
//                      comparisonFun = "==",
//                      limitConstraints = 12
//                  ) -->> Currently no function for this case. We do the same
//                   thing as we do for the CoarseGrained case
//
// PrmDstPartNoZ  : permuteCount(
//                      1:20, 4,
//                      constraintFun = "sum",
//                      comparisonFun = "==",
//                      limitConstraints = 20
//                  ) -->> CountCompsDistinctLen(20, 4)
//
//                  ## When zero is involved, we call the same compiled
//                  ## function just translated by the width. Note in this
//                  ## case, we are passing includeZero = FALSE b/c there
//                  ## is only one zero and hence isomorphic to 1:24
// PrmDstPrtOneZ  : permuteCount(
//                      0:20, 4,
//                      constraintFun = "sum",
//                      comparisonFun = "==",
//                      limitConstraints = 20
//                  ) -->> CountCompsDistinctLen(24, 4)
//
// PrmDstPartMZ   : permuteCount(
//                      0:20, 4,
//                      freqs = c(2, rep(1, 20)),
//                      constraintFun = "sum",
//                      comparisonFun = "==",
//                      limitConstraints = 20
//                  ) -->> CountCompsDistinctMZWeak(20, 4, 20, 2)
//
// double CountPartsPermDistinctCap(const std::vector<int> &z, int cap,
//                                  int tar, int width, bool includeZero)
//
// For the cases below, we point out that the sum(z) will not be equal to the
// target as z in these cases is being utilized as the index vector. These
// cases are under the "general" case and require a vector, v, that is used
// for output. E.g. in all of the generating functions, you will see v[z[i]].
//
// PrmDstPrtCap   : permuteCount(
//                      55, 4,
//                      constraintFun = "sum",
//                      comparisonFun = "==",
//                      limitConstraints = 80
//                  ) -->> CountPartsPermDistinctCap(80, 4, 55)
//
//                  ## When zero is involved, we call the same compiled
//                  ## function just translated by the width. Note in this
//                  ## case, we are passing includeZero = FALSE b/c there
//                  ## is only one zero and hence isomorphic to 1:84
//                  permuteCount(
//                      0:55, 4,
//                      constraintFun = "sum",
//                      comparisonFun = "==",
//                      limitConstraints = 80
//                  ) -->> CountPartsPermDistinctCap(84, 4, 56)
//
// PrmDstPrtCapMZ : permuteCount(
//                      0:55, 4,
//                      freqs = c(2, rep(1, 55)),
//                      constraintFun = "sum",
//                      comparisonFun = "==",
//                      limitConstraints = 80
//                  ) -->> CountPartsPermDistinctCap(80, 4, 55, 2)
//
// NotPartition   :
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
// The cases below are technically compositions, however we don't have a
// composition algorithm for them. We do have a partition algorithm available
// for these cases. Remember, in order for us to have a partition/composition
// algorithm, it must satisfy the requirement of producing the next
// lexicographical result. These cases rely on generating the next partition
// and subsequently generating all possible permutations of such partition.
// This will not produce the results in lexicographical order, but it will
// produce all results.
//
// Note, if zero is included, it will be considered when generating perms,
// thus all cases below will produce weak compositions.
//
// Also note that all of the cases below will stem from permuteCount/General
//
// PrmRepPartNoZ  : Permutations of partitions with repetition & no zeros
// PrmRepPart     : Permutations of partitions with repetition
// PrmRepCapped   : Perms of partitions with repetition & restricted parts
// PrmDstPartNoZ  : Permutations of partitions with distinct parts & no zeros
// PrmDstPrtOneZ  : Permutations of partitions with distinct parts & one zero
// PrmDstPartMZ   : Permutations of partitions with distinct parts & more
//                   than one zero
//
// PrmDstPrtCap   : Permutations of partitions with distinct & restricted parts
// PrmDstPrtCapMZ : Permutations of partitions with distinct & restricted parts
//                   and more than one zero
//
// PrmMultiset    : Permutations of partitions of non-trivial multisets.
//                   Non-trivial here means elements other than 0 have
//                   multiplicity > 1
//
// NotMapped      : These are partitions, however they are not mapped
// NoSolution     : This passes the CheckPartition function, but there is no
//                   solution given the target, width, or constraints
// NotPartition   : Does not pass the CheckPartition function
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
    PrmRepPartNoZ  = 20,
    PrmRepPart     = 21,
    PrmRepCapped   = 22,
    PrmDstPartNoZ  = 23,
    PrmDstPrtOneZ  = 24,
    PrmDstPartMZ   = 25,
    PrmDstPrtCap   = 26,
    PrmDstPrtCapMZ = 27,
    PrmMultiset    = 28,
    NotMapped      = 29,
    NoSolution     = 30,
    NotPartition   = 31,
    NumTypes       = 32
};

constexpr const char* PTypeNames[] = {
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
    "PrmRepPartNoZ",
    "PrmRepPart",
    "PrmRepCapped",
    "PrmDstPartNoZ",
    "PrmDstPrtOneZ",
    "PrmDstPartMZ",
    "PrmDstPrtCap",
    "PrmDstPrtCapMZ",
    "PrmMultiset",
    "NotMapped",
    "NoSolution",
    "NotPartition",
    "NumTypes"
};

const std::array<PartitionType, 5> NoCountAlgoPTypeArr{{
    PartitionType::NotMapped, PartitionType::PrmRepCapped,
    PartitionType::NotPartition, PartitionType::NoSolution,
    PartitionType::CoarseGrained
}};

const std::array<PartitionType, 6> CappedPTypeArr{{
    PartitionType::RepCapped, PartitionType::DstctCapped,
    PartitionType::DstctCappedMZ, PartitionType::PrmRepCapped,
    PartitionType::PrmDstPrtCap, PartitionType::PrmDstPrtCapMZ
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
    bool isPerm = false;      // This is to distinguish between true
                              // integer compositions and permutations of
                              // integer partitions. The latter occurs when we
                              // have an algorithm for generating a particular
                              // type of partition but no known algorithm for
                              // generating compositions.
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
