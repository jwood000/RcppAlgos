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
// RepStdAll       : CountPartRep(20)
// RepNoZero       : CountPartRepLen(20, 5)
// RepShort        : CountPartRepLen(23, 3)
// RepCapped       : CountPartRepLenCap(14, 3, 10) N.B. Get first part:
//                    (3, 5, 12); Map to match(c(3, 5, 12), 3:12) -->>
//                    (1, 3, 10) -->> sum(c(1, 3, 10)) = 14
//
// DstctStdAll     : CountPartDistinct(20)
// DstctMultiZero  : CountPartsDistinctMultiZero(c(0, 0, 1, 2, 17), 20, 5)
// DstctOneZero    : CountPartDistinctLen(25, 5) N.B. Add 1 to each element to
//                    obtain new target = 25
//
// DstctNoZero     : CountPartDistinctLen(20, 5)
// DstctCapped     : CountPartDistinctLenCap(20, 4, 9)
// DstctCappedMZ   : CountPartsDistinctMultiZeroCap(c(0, 0, 9, 11), 20, 4, 11)
// LengthOne       : 1 or 0
// Multiset        : CountPartsMultiset(rep(1:3, 5), c(1, 2, 2, 15))
//
// CoarseGrained   : Currently, we utilize std::vector and push_back until we
//                    reach a terminating situation.
//
// CompRepNoZero   : compositionsCount(20, 5, TRUE) -->> CountCompsRepLen(20, 5)
// CompRepWeak     : compositionsCount(0:20, 5, TRUE, weak = TRUE) -->>
//                    CountCompsRepLen(25, 5) We add the length to the target
//
// CmpRpZroNotWk   : compositionsCount(0:20, 5, TRUE) -- >>
//                    CountCompsRepZNotWk(20, 5)
//
// CmpDstctNoZero  : compositionsCount(20, 5) -->> CountCompsDistinctLen(20, 5)
// CmpDstctZNotWk  : compositionsCount(0:20, 5, freqs = c(3, rep(1, 20))) -->>
//                    CountCompsDistinctMultiZero(20, 5, 0, 2)
//
// CmpDstctWeak    : compositionsCount(0:20, 5, weak = TRUE) -->>
//                    CountCompsDistinctLen(25, 5)
//
// CmpDstctMZWeak  : compositionsCount(0:20, 5, freqs = c(3, rep(1, 20)),
//                                    weak = TRUE) -->>
//                    CountCompsDistinctMZWeak(20, 5, 0, 2)
//
// CmpDstctCapped  : compositionsCount(10, 4, target = 25) -->>
//                    CountCompDistLenRstrctd(25, 4, {1, 2, ..., 10})
//
// CmpDstCapWeak   : compositionsCount(0:10, 4, target = 25, weak = TRUE) -->>
//                    CountCompDistLenRstrctd(25, 4, {1, 2, ..., 10})
//
// CmpDstCapMZNotWk: compositionsCount(0:10, 4, target = 25) -->>
//                    CountCompsDistinctRstrctdMZ(25, 4, {1, 2, ..., 10}, 3)
//
//                   compositionsCount(0:13, 4, freqs = c(2, rep(1, 13)),
//                                    target = 25) -->>
//                    CountCompsDistinctRstrctdMZ(25, 4, {1, 2, ..., 13}, 2)
//
// CmpDstCapMZWeak : compositionsCount(0:20, 5, freqs = c(4, rep(1, 20)),
//                                     target = 35, weak = TRUE)
//                    ## N.B. Even though we could have 4 zeros, it is
//                    ## automatically determined that the shortest len is 2.
//                    CountPartsPermDistinctRstrctdMZ(35, 5, {1..20}, 2);
//
// CompMultiset    : CountPartsMultiset(rep(1:3, 5), c(1, 2, 2, 15), true)
//
// PrmRepPartNoZ   : permuteCount(
//                       1:20, 5, TRUE,
//                       constraintFun = "sum",
//                       comparisonFun = "==",
//                       limitConstraints = 20
//                   ) -->> CountCompsRepLen(20, 5)
//
//                   ## When zero is involved, we call the same compiled
//                   ## function just translated by the width
// PrmRepPart      : permuteCount(
//                       0:20, 5, TRUE,
//                       constraintFun = "sum",
//                       comparisonFun = "==",
//                       limitConstraints = 20
//                   ) -->> CountCompsRepLen(25, 5)
//
// PrmRepCapped    : permuteCount(
//                       1:7, 5, TRUE,
//                       constraintFun = "sum",
//                       comparisonFun = "==",
//                       limitConstraints = 12
//                   ) -->> Currently no function for this case. We do the same
//                    thing as we do for the CoarseGrained case
//
// PrmDstPartNoZ   : permuteCount(
//                       1:20, 4,
//                       constraintFun = "sum",
//                       comparisonFun = "==",
//                       limitConstraints = 20
//                   ) -->> CountCompsDistinctLen(20, 4)
//
//                   ## When zero is involved, we call the same compiled
//                   ## function just translated by the width. Note in this
//                   ## case, we are passing includeZero = FALSE b/c there
//                   ## is only one zero and hence isomorphic to 1:24
// PrmDstPrtOneZ   : permuteCount(
//                       0:20, 4,
//                       constraintFun = "sum",
//                       comparisonFun = "==",
//                       limitConstraints = 20
//                   ) -->> CountCompsDistinctLen(24, 4)
//
// PrmDstPartMZ    : permuteCount(
//                       0:20, 4,
//                       freqs = c(2, rep(1, 20)),
//                       constraintFun = "sum",
//                       comparisonFun = "==",
//                       limitConstraints = 20
//                   ) -->> CountCompsDistinctMZWeak(20, 4, 20, 2)
//
// double CountPartsPermDistinctCap(const std::vector<int> &z, int cap,
//                                  int tar, int width, bool includeZero)
//
// For the cases below, we point out that the sum(z) will not be equal to the
// target as z in these cases is being utilized as the index vector. These
// cases are under the "general" case and require a vector, v, that is used
// for output. E.g. in all of the generating functions, you will see v[z[i]].
//
// PrmDstPrtCap    : permuteCount(
//                       55, 4,
//                       constraintFun = "sum",
//                       comparisonFun = "==",
//                       limitConstraints = 80
//                   ) -->> CountPartsPermDistinctCap(80, 4, 55)
//
//                   ## When zero is involved, we call the same compiled
//                   ## function just translated by the width. Note in this
//                   ## case, we are passing includeZero = FALSE b/c there
//                   ## is only one zero and hence isomorphic to 1:84
//                   permuteCount(
//                       0:55, 4,
//                       constraintFun = "sum",
//                       comparisonFun = "==",
//                       limitConstraints = 80
//                   ) -->> CountPartsPermDistinctCap(84, 4, 56)
//
// PrmDstPrtCapMZ  : permuteCount(
//                       0:55, 4,
//                       freqs = c(2, rep(1, 55)),
//                       constraintFun = "sum",
//                       comparisonFun = "==",
//                       limitConstraints = 80
//                   ) -->> CountPartsPermDistinctCap(80, 4, 55, 2)
//
// PrmMultiset     : permuteCount(
//                        1:20, 5, freqs = rep(1:4, 5), limitConstraints = 25,
//                        constraintFun = "sum", comparisonFun = "=="
//                   ) -->> CountPartsMultiset(
//                              rep(1:4, 5), {1, 2, 2, 3, 17}, true, true
//                          )
//
// NotMapped       :
// NoSolution      :
// NotPartition    :
//
// ****************************************************************************
//
//
// ************************** Definitions w/ Examples *************************
//
// Notes:
// * startZ is the canonical "first" index/result vector in the mapped/standard
//   problem space used by the core next-lex algorithms.
// * "Capped" means parts are restricted to a finite window of v (i.e. cap).
// * "MZ" (MultiZero) refers to cases where 0 may appear multiple times due to
//   freqs[0] (or equivalent mapping), and startZ may or may not maximize 0s.
//
// RepStdAll       : Get all partitions with repetition (0 allowed if present).
//                   E.g. tar = 20; startZ = c(0, 0, 0, 0, 20):
//                   CountPartRep(20)
//
// RepNoZero       : Partitions with repetition excluding 0.
//                   E.g. tar = 20; m = 5; startZ = c(1, 1, 1, 1, 16):
//                   CountPartRepLen(20, 5)
//
// RepShort        : Repetition case where width is not maximized (m fixed,
//                   m < tar in the include-zero design).
//                   E.g. tar = 20; m = 3; startZ = c(0, 0, 20)
//
// RepCapped       : Repetition partitions with restricted parts.
//                   E.g. tar = 20, m = 3, v = 3:12:
//                   (mapped tar = 14 from 0:9) startZ ~ c(0, 0, 14)
//
// DstctStdAll     : Get all distinct partitions (0 may repeat via freqs[0]).
//                   E.g. tar = 20; startZ = c(0, 0, 0, 0, 20)
//
// DstctMultiZero  : Distinct parts where multiple zeros are possible; startZ
//                   does not necessarily maximize the number of zeros.
//                   E.g. tar = 20; startZ = c(0, 0, 1, 2, 17)
//
// DstctOneZero    : Distinct parts where exactly one zero is possible.
//                   Often seen when isMult = FALSE but 0 is included.
//                   E.g. tar = 20; startZ = c(0, 1, 2, 3, 14)
//
// DstctNoZero     : Distinct partitions excluding 0.
//                   E.g. tar = 20; startZ = c(1, 2, 3, 4, 10)
//
// DstctCapped     : Distinct partitions with restricted parts (cap/window).
//                   E.g. tar = 20, m = 4, v = 1:9 gives startZ = c(1, 2, 8, 9)
//
// DstctCappedMZ   : Distinct + capped + multi-zero (freqs[0] > 1).
//                   E.g. tar = 20, m = 4, v = 0:11, freqs = c(2, rep(1, 11))
//                   gives startZ = c(0, 0, 9, 11)
//
// LengthOne       : Any partition/composition when m = 1.
//
// Multiset        : Partitions of non-trivial multisets. Non-trivial here means
//                   elements other than 0 have multiplicity > 1.
//
// CoarseGrained   : Partition-esque constraints that pass CheckPartition but
//                   do not admit a dedicated next-lex partition algorithm.
//                   This is equal to ConstraintType::PartitionEsque.
//
// CompRepNoZero   : Standard compositions with repetition and no zeros.
//                   E.g. tar = 20, m = 5; startZ = c(1, 1, 1, 1, 16)
//
// CompRepWeak     : Repetition compositions where zeros are allowed (weak).
//                   E.g. tar = 20, m = 5; startZ = c(0, 0, 0, 0, 20)
//
// CmpRpZroNotWk   : Compositions where 0 is in v, but we only want permutations
//                   of non-zero values (non-weak output). Internally behaves
//                   like repetition comps with a "zero slot" used for mapping.
//                   E.g. tar = 20, m = 5; startZ = c(0, 0, 0, 0, 20)
//
// CmpDstctNoZero  : Standard compositions with distinct parts and no zeros.
//                   E.g. tar = 20; m = 5; startZ = c(1, 2, 3, 4, 10)
//
// CmpDstctZNotWk  : Distinct compositions where 0 exists/may appear, but only
//                   non-zero values participate in next-iteration mechanics.
//                   (Conceptually: distinct, non-weak, with a mapped zero slot.)
//                   E.g. tar = 20; m = 5; startZ = c(0, 1, 2, 3, 14)
//
// CmpDstctWeak    : Distinct weak compositions where a single zero is allowed.
//                   E.g. tar = 20; m = 5; startZ = c(0, 1, 2, 3, 14)
//
// CmpDstctMZWeak  : Distinct weak compositions where multiple zeros are allowed
//                   (freqs[0] > 1 / multi-zero mapping).
//                   E.g. tar = 20; m = 5; startZ = c(0, 0, 1, 2, 17)
//
// CmpDstctCapped  : Distinct compositions with restricted parts (cap/window).
//                   E.g. tar = 20, m = 4, v = 1:9 gives startZ = c(1, 2, 8, 9)
//
// CmpDstCapWeak   : Distinct capped weak compositions (0 allowed, capped set).
//                   E.g. tar = 20, m = 4, v = 0:9 gives startZ = c(0, 1, 8, 11)
//                   (mapped startZ shown; exact starter depends on cap/target)
//
// CmpDstCapMZNotWk: Distinct capped + multi-zero, non-weak iteration rules.
//                   E.g. tar = 20, m = 4, v = 0:11, freqs = c(2, rep(1, 11))
//                   gives startZ = c(0, 0, 9, 11)
//
// CmpDstCapMZWeak : Same as above, but weak (zeros allowed in results).
//
// CompMultiset    : Compositions of non-trivial multisets. Non-trivial here
//                   means elements other than 0 have multiplicity > 1.
//
// The cases below are technically compositions, however we do not have a
// dedicated next-lex composition algorithm for them. We instead generate the
// next partition and then enumerate its permutations. This produces all
// results but not in lexicographical order.
//
// Note: If zero is included, it is considered in permutation generation, so
// these produce weak compositions when 0 is present.
//
// PrmRepPartNoZ   : Permutations of repetition partitions with no zeros.
// PrmRepPart      : Permutations of repetition partitions (0 may appear).
// PrmRepCapped    : Permutations of repetition partitions with restricted parts.
// PrmDstPartNoZ   : Permutations of distinct partitions with no zeros.
// PrmDstPrtOneZ   : Permutations of distinct partitions with exactly one zero.
// PrmDstPartMZ    : Permutations of distinct partitions with multiple zeros.
// PrmDstPrtCap    : Permutations of distinct capped partitions.
// PrmDstPrtCapMZ  : Permutations of distinct capped partitions with multi-zero.
// PrmMultiset     : Permutations of partitions of non-trivial multisets.
//
// NotMapped       : Partition-like input, but mapping heuristics did not
//                   identify an isomorphic standard/capped case.
// NoSolution      : Passes CheckPartition, but no solution exists for the
//                   given target, width, and/or constraints.
// NotPartition    : Does not pass CheckPartition.
//
// ****************************************************************************

enum class PartitionType {
    RepStdAll        = 0,
    RepNoZero        = 1,
    RepShort         = 2,
    RepCapped        = 3,
    DstctStdAll      = 4,
    DstctMultiZero   = 5,
    DstctOneZero     = 6,
    DstctNoZero      = 7,
    DstctCapped      = 8,
    DstctCappedMZ    = 9,
    LengthOne        = 10,
    Multiset         = 11,
    CoarseGrained    = 12,
    CompRepNoZero    = 13,
    CompRepWeak      = 14,
    CmpRpZroNotWk    = 15,
    CmpDstctNoZero   = 16,
    CmpDstctZNotWk   = 17,
    CmpDstctWeak     = 18,
    CmpDstctMZWeak   = 19,
    CmpDstctCapped   = 20,
    CmpDstCapWeak    = 21,
    CmpDstCapMZNotWk = 22,
    CmpDstCapMZWeak  = 23,
    CompMultiset     = 24,
    PrmRepPartNoZ    = 25,
    PrmRepPart       = 26,
    PrmRepCapped     = 27,
    PrmDstPartNoZ    = 28,
    PrmDstPrtOneZ    = 29,
    PrmDstPartMZ     = 30,
    PrmDstPrtCap     = 31,
    PrmDstPrtCapMZ   = 32,
    PrmMultiset      = 33,
    NotMapped        = 34,
    NoSolution       = 35,
    NotPartition     = 36,
    NumTypes         = 37
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
    "CmpDstctWeak",
    "CmpDstctMZWeak",
    "CmpDstctCapped",
    "CmpDstCapWeak",
    "CmpDstCapMZNotWk",
    "CmpDstCapMZWeak",
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

const std::array<PartitionType, 10> CappedPTypeArr{{
    PartitionType::RepCapped, PartitionType::DstctCapped,
    PartitionType::DstctCappedMZ, PartitionType::PrmRepCapped,
    PartitionType::PrmDstPrtCap, PartitionType::PrmDstPrtCapMZ,
    PartitionType::CmpDstctCapped, PartitionType::CmpDstCapWeak,
    PartitionType::CmpDstCapMZWeak, PartitionType::CmpDstCapMZNotWk
}};

const std::array<PartitionType, 8> CmpDstPTypeArr{{
    PartitionType::CmpDstctWeak, PartitionType::CmpDstCapWeak,
    PartitionType::CmpDstctMZWeak, PartitionType::CmpDstctNoZero,
    PartitionType::CmpDstctZNotWk, PartitionType::CmpDstctCapped,
    PartitionType::CmpDstCapMZWeak, PartitionType::CmpDstCapMZNotWk
}};

struct PartDesign {
    // Maximum number of zeros permitted in generated results.
    // Needed to reconstruct the complement vector when iteration
    // starts from an arbitrary result that may currently contain no zeros.
    int maxZeros = 0;
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
