#pragma once

#include <cstdint>
#include <vector>
#include <gmpxx.h>
#include <array>

// ************************** Count Function Examples *************************
//
// The examples below show the counting routine associated with each
// PartitionType. Some cases require a mapped target and/or mapped set of
// allowed values before calling the underlying count routine.
//
// RepStdAll       : CountPartRep(20)
//
// RepNoZero       : CountPartRepLen(20, 5)
//
// RepShort        : CountPartRepLen(23, 3)
//
// RepCapped       : CountPartRepLenCap(14, 3, 10)
//                   N.B. First partition: (3, 5, 12)
//                   Map with match(c(3, 5, 12), 3:12) -> (1, 3, 10)
//                   giving mapped target 14.
//
// DistinctStdAll  : CountPartDistinct(20)
//
// DistinctMZ      : CountPartsDistinctMZ(
//                       c(0, 0, 1, 2, 17), 20, 5
//                   )
//
// DistinctOneZero : CountPartDistinctLen(25, 5)
//                   N.B. Add 1 to each part, giving mapped target 25.
//
// DistinctNoZero  : CountPartDistinctLen(20, 5)
//
// DistinctCapped  : CountPartDistinctLenCap(20, 4, 9)
//
// DistinctCappedMZ
//                 : CountPartsDistinctMZCap(
//                       c(0, 0, 9, 11), 20, 4, 11
//                   )
//
// LengthOne       : 1 or 0
//
// Multiset        : CountPartsMultiset(20, 4, rep(1:3, 5))
//
// CoarseGrained   : No dedicated count routine. Generation uses a
//                   std::vector and push_back until the terminating
//                   condition is reached.
//
// CompRepNoZero   : compositionsCount(20, 5, TRUE)
//                   -> CountCompsRepLen(20, 5)
//
// CompRepWeak     : compositionsCount(0:20, 5, TRUE, weak = TRUE)
//                   -> CountCompsRepLen(25, 5)
//                   N.B. Add the width to the target.
//
// CompRepWeakCap  : compositionsCount(
//                       0:20, 5, TRUE, weak = TRUE, target = 40
//                   )
//                   -> CountCompsRepLenCap(45, 5, 1:20)
//                   N.B. Add the width to the target.
//
// CompRepCapped   : compositionsCount(
//                       3, 6, repetition = TRUE, target = 10
//                   )
//                   -> CountCompsRepLenCap(10, 6, 1:3)
//
// CompRepCapZero  : compositionsCount(
//                       0:3, 6, repetition = TRUE, target = 10
//                   )
//                   -> CountCompsRepCapZero(10, 6, 1:3)
//
// CompRepZero     : compositionsCount(0:20, 5, TRUE)
//                   -> CountCompsRepZero(20, 5)
//
// CompDistinctNoZero
//                 : compositionsCount(20, 5)
//                   -> CountCompsDistinctLen(20, 5)
//
// CompDistinctZero
//                 : compositionsCount(
//                       0:20, 5, freqs = c(3, rep(1, 20))
//                   )
//                   -> CountCompsDistinctMZ(20, 5, 0, 2)
//
// CompDistinctWeak
//                 : compositionsCount(0:20, 5, weak = TRUE)
//                   -> CountCompsDistinctLen(25, 5)
//
// CompDistinctMZWeak
//                 : compositionsCount(
//                       0:20, 5,
//                       freqs = c(3, rep(1, 20)),
//                       weak = TRUE
//                   )
//                   -> CountCompsDistinctMZWeak(20, 5, 0, 2)
//
// CompDistinctCapped
//                 : compositionsCount(10, 4, target = 25)
//                   -> CountCompsDistLenRstrctd(
//                          25, 4, {1, 2, ..., 10}
//                      )
//
// CompDistinctCapWeak
//                 : compositionsCount(
//                       0:10, 4, target = 25, weak = TRUE
//                   )
//                   -> CountCompsDistLenRstrctd(
//                          25, 4, {1, 2, ..., 10}
//                      )
//
// CompDistinctCapMZ
//                 : compositionsCount(0:10, 4, target = 25)
//                   -> CountCompsDistinctRstrctdMZ(
//                          25, 4, {1, 2, ..., 10}, 3
//                      )
//
//                   compositionsCount(
//                       0:13, 4,
//                       freqs = c(2, rep(1, 13)),
//                       target = 25
//                   )
//                   -> CountCompsDistinctRstrctdMZ(
//                          25, 4, {1, 2, ..., 13}, 2
//                      )
//
// CompDistinctCapMZWeak
//                 : compositionsCount(
//                       0:20, 5,
//                       freqs = c(4, rep(1, 20)),
//                       target = 35,
//                       weak = TRUE
//                   )
//                   N.B. Although four zeros are available, the shortest
//                   feasible positive length is determined automatically.
//                   -> CountPartsPermDistinctRstrctdMZ(
//                          35, 5, {1, ..., 20}, 2
//                      )
//
// CompMultiset    : CountCompsMultiset(20, 4, rep(1:3, 5))
//
// CompMultisetZero
//                 : Non-weak multiset composition count with at least 1 mapped
//                   zeros. The count spans the feasible positive composition
//                   lengths represented by those zero-padding slots.
//
// CompMultisetWeak
//                 : Weak multiset composition count with at least one zero.
//                   Zero participates as an actual composition part.
//
// PrmRepPartNoZ   : permuteCount(
//                       1:20, 5, TRUE,
//                       constraintFun = "sum",
//                       comparisonFun = "==",
//                       limitConstraints = 20
//                   )
//                   -> CountCompsRepLen(20, 5)
//
//                   When zero is involved, the same compiled routine is used
//                   after translating the target by the width.
//
// PrmRepPart      : permuteCount(
//                       0:20, 5, TRUE,
//                       constraintFun = "sum",
//                       comparisonFun = "==",
//                       limitConstraints = 20
//                   )
//                   -> CountCompsRepLen(25, 5)
//
// PrmRepCapped    : permuteCount(
//                       1:7, 5, TRUE,
//                       constraintFun = "sum",
//                       comparisonFun = "==",
//                       limitConstraints = 12
//                   )
//                   No dedicated count routine currently exists. This case
//                   uses the same dynamic-generation strategy as
//                   CoarseGrained.
//
// PrmDstPartNoZ   : permuteCount(
//                       1:20, 4,
//                       constraintFun = "sum",
//                       comparisonFun = "==",
//                       limitConstraints = 20
//                   )
//                   -> CountCompsDistinctLen(20, 4)
//
//                   When exactly one zero is involved, the problem is
//                   translated to an isomorphic positive-valued problem.
//
// PrmDstPrtOneZ   : permuteCount(
//                       0:20, 4,
//                       constraintFun = "sum",
//                       comparisonFun = "==",
//                       limitConstraints = 20
//                   )
//                   -> CountCompsDistinctLen(24, 4)
//
// PrmDstPartMZ    : permuteCount(
//                       0:20, 4,
//                       freqs = c(2, rep(1, 20)),
//                       constraintFun = "sum",
//                       comparisonFun = "==",
//                       limitConstraints = 20
//                   )
//                   -> CountCompsDistinctMZWeak(20, 4, 20, 2)
//
// For the capped permutation cases below, sum(z) need not equal the original
// target because z is an index vector into v. These are general/mapped cases;
// generation obtains output values through expressions such as v[z[i]].
//
// PrmDstPrtCap    : permuteCount(
//                       55, 4,
//                       constraintFun = "sum",
//                       comparisonFun = "==",
//                       limitConstraints = 80
//                   )
//                   -> CountPartsPermDistinctCap(80, 4, 55)
//
//                   With exactly one zero, the problem is translated to an
//                   isomorphic positive-valued problem:
//
//                   permuteCount(
//                       0:55, 4,
//                       constraintFun = "sum",
//                       comparisonFun = "==",
//                       limitConstraints = 80
//                   )
//                   -> CountPartsPermDistinctCap(84, 4, 56)
//
// PrmDstPrtCapMZ  : permuteCount(
//                       0:55, 4,
//                       freqs = c(2, rep(1, 55)),
//                       constraintFun = "sum",
//                       comparisonFun = "==",
//                       limitConstraints = 80
//                   )
//                   -> CountPartsPermDistinctCap(80, 4, 55, 2)
//
// PrmMultiset     : permuteCount(
//                       1:20, 5,
//                       freqs = rep(1:4, 5),
//                       limitConstraints = 25,
//                       constraintFun = "sum",
//                       comparisonFun = "=="
//                   )
//                   -> CountPartsMultiset(
//                          rep(1:4, 5),
//                          {1, 2, 2, 3, 17},
//                          true,
//                          true
//                      )
//
// NotMapped       : No dedicated partition count routine.
//
// NoSolution      : Count is zero.
//
// NotPartition    : Uses the general constraint machinery rather than the
//                   partition-specific count routines.
//
// ****************************************************************************
//
// ************************** Definitions w/ Examples *************************
//
// Notes:
//
// * startZ is the canonical first index/result vector in the mapped or standard
//   problem space used by the core next-lex algorithms.
//
// * "Capped" means parts are restricted to a finite window of v (i.e. cap).
//
// * "MZ" (MZ) means multiple zeros may be available through freqs[0]
//   or an equivalent mapped representation. MZ does not imply weakness.
//
// * Composition types are non-weak unless "Weak" is explicitly present in the
//   PartitionType name.
//
// * For non-weak compositions, zeros appearing in startZ are internal padding
//   used to represent compositions having fewer than width positive parts.
//   Those zeros are not returned as composition parts.
//
// * For weak compositions, zero is an actual part. It participates in the
//   result, counting, generation, ranking, and ordering semantics.
//
// RepStdAll       : All repetition partitions in the standard design. Zero may
//                   appear when included by the design.
//                   E.g. tar = 20;
//                   startZ = c(0, 0, 0, 0, 20)
//
// RepNoZero       : Fixed-width repetition partitions excluding zero.
//                   E.g. tar = 20; m = 5;
//                   startZ = c(1, 1, 1, 1, 16)
//
// RepShort        : Fixed-width repetition partitions where width is smaller
//                   than the maximal include-zero design width.
//                   E.g. tar = 20; m = 3;
//                   startZ = c(0, 0, 20)
//
// RepCapped       : Repetition partitions with parts restricted by cap/window.
//                   E.g. tar = 20; m = 3; v = 3:12;
//                   mapped tar = 14;
//                   startZ ~ c(0, 0, 14)
//
// DistinctStdAll  : Partitions whose non-zero parts are distinct, with zero
//                   allowed to repeat through freqs[0] when applicable.
//                   E.g. tar = 20;
//                   startZ = c(0, 0, 0, 0, 20)
//
// DistinctMZ      : Partitions with distinct non-zero parts and multiple zeros
//                   available. startZ need not maximize the number of zeros.
//                   E.g. tar = 20;
//                   startZ = c(0, 0, 1, 2, 17)
//
// DistinctOneZero : Distinct partitions where at most one zero is available.
//                   Often encountered when isMult = FALSE and zero is present.
//                   E.g. tar = 20;
//                   startZ = c(0, 1, 2, 3, 14)
//
// DistinctNoZero  : Distinct partitions excluding zero.
//                   E.g. tar = 20;
//                   startZ = c(1, 2, 3, 4, 10)
//
// DistinctCapped  : Distinct partitions with parts restricted by cap/window.
//                   E.g. tar = 20; m = 4; v = 1:9;
//                   startZ = c(1, 2, 8, 9)
//
// DistinctCappedMZ
//                 : Capped partitions with distinct non-zero parts and
//                   multiple zeros available.
//                   E.g. tar = 20; m = 4; v = 0:11;
//                   freqs = c(2, rep(1, 11));
//                   startZ = c(0, 0, 9, 11)
//
// LengthOne       : Any partition/composition problem whose width is one.
//
// Multiset        : Partitions of a non-trivial multiset without a multi-zero
//                   mapped state. "Non-trivial" means at least one non-zero
//                   value has multiplicity greater than one.
//
// CoarseGrained   : Partition-like constraints that pass CheckPartition but do
//                   not admit a dedicated next-lex partition algorithm.
//                   Corresponds to ConstraintType::PartitionEsque.
//
// CompRepNoZero   : Standard fixed-width repetition compositions excluding
//                   zero.
//                   E.g. tar = 20; m = 5;
//                   startZ = c(1, 1, 1, 1, 16)
//
// CompRepWeak     : Weak repetition compositions. Zero is an actual part.
//                   E.g. tar = 20; m = 5;
//                   startZ = c(0, 0, 0, 0, 20)
//
// CompRepWeakCap  : Capped weak repetition compositions. Zero is an actual
//                   composition part.
//                   E.g. tar = 40; m = 5; cap = 20;
//                   startZ = c(0, 0, 0, 20, 20)
//
// CompRepCapped   : Capped repetition compositions excluding zero-padding.
//                   E.g. tar = 10; m = 5; cap = 3;
//                   startZ = c(1, 1, 2, 3, 3)
//
// CompRepCapZero  : Capped non-weak repetition compositions with zero present
//                   in the mapped problem. Zero acts only as internal padding
//                   for shorter positive compositions and is not returned as
//                   a composition part.
//
// CompRepZero     : Non-weak repetition compositions with zero present in the
//                   mapped problem. Zero acts only as internal padding for
//                   shorter positive compositions and is not returned as a
//                   composition part.
//                   E.g. tar = 20; m = 5;
//                   startZ = c(0, 0, 0, 0, 20)
//
// CompDistinctNoZero
//                 : Standard compositions with distinct parts and no zero.
//                   E.g. tar = 20; m = 5;
//                   startZ = c(1, 2, 3, 4, 10)
//
// CompDistinctZero
//                 : Non-weak compositions with distinct positive parts and a
//                   mapped zero slot. Zero is internal padding and does not
//                   participate in the returned composition.
//                   E.g. tar = 20; m = 5;
//                   startZ = c(0, 1, 2, 3, 14)
//
// CompDistinctWeak
//                 : Weak compositions with distinct positive parts and at most
//                   one zero. Zero is an actual composition part.
//                   E.g. tar = 20; m = 5;
//                   startZ = c(0, 1, 2, 3, 14)
//
// CompDistinctMZWeak
//                 : Weak compositions with distinct positive parts and
//                   multiple zeros available. Zeros are actual composition
//                   parts.
//                   E.g. tar = 20; m = 5;
//                   startZ = c(0, 0, 1, 2, 17)
//
// CompDistinctCapped
//                 : Capped compositions with distinct positive parts and no
//                   zero-padding.
//                   E.g. tar = 20; m = 4; v = 1:9;
//                   startZ = c(1, 2, 8, 9)
//
// CompDistinctCapWeak
//                 : Capped weak compositions with distinct positive parts and
//                   zero available as an actual composition part.
//
// CompDistinctCapMZ
//                 : Capped non-weak compositions with distinct positive parts
//                   and multiple mapped zeros available. Zeros act only as
//                   internal padding and do not participate in the returned
//                   composition.
//                   E.g. tar = 20; m = 4; v = 0:11;
//                   freqs = c(2, rep(1, 11));
//                   startZ = c(0, 0, 9, 11)
//
// CompDistinctCapMZWeak
//                 : Capped weak compositions with distinct positive parts and
//                   multiple zeros available. Zeros are actual composition
//                   parts.
//
// CompMultiset    : Compositions of a non-trivial multiset without a
//                   multi-zero mapped state. At least one non-zero value has
//                   multiplicity greater than one.
//
// CompMultisetZero
//                 : Non-weak compositions of a non-trivial multiset with
//                   multiple zeros available in the mapped/design state.
//                   Zeros are internal padding representing compositions with
//                   fewer than width positive parts and are not returned.
//
// CompMultisetWeak
//                 : Weak compositions of a non-trivial multiset with multiple
//                   zeros available. Zero is an actual composition part and
//                   participates in counting, generation, ranking, and
//                   ordering.
//
// The types below are composition-like results produced through partition
// generation followed by permutation generation. We do not have dedicated
// next-lex composition algorithms for these cases.
//
// The resulting order is therefore not guaranteed to be lexicographical.
//
// If zero is present in the generated partition, it is treated as an ordinary
// part during permutation generation. Consequently, there is no separate
// weak/non-weak distinction for these permutation-based types: a result
// containing zero is inherently weak.
//
// PrmRepPartNoZ   : Permutations of repetition partitions containing no zero.
//
// PrmRepPart      : Permutations of repetition partitions where zero may be
//                   present.
//
// PrmRepCapped    : Permutations of capped repetition partitions.
//
// PrmDstPartNoZ   : Permutations of distinct partitions containing no zero.
//
// PrmDstPrtOneZ   : Permutations of distinct partitions containing at most
//                   one zero.
//
// PrmDstPartMZ    : Permutations of partitions with distinct positive parts
//                   and multiple zeros.
//
// PrmDstPrtCap    : Permutations of capped distinct partitions.
//
// PrmDstPrtCapMZ  : Permutations of capped partitions with distinct positive
//                   parts and multiple zeros.
//
// PrmMultiset     : Permutations of partitions of non-trivial multisets. Zero
//                   is always treated as an ordinary part when present. There
//                   is no non-weak permutation-multiset variant.
//
// NotMapped       : Partition-like input for which the mapping heuristics did
//                   not identify an isomorphic standard/capped case.
//
// NoSolution      : Input passes CheckPartition, but no solution exists for
//                   the target, width, and/or supplied constraints.
//
// NotPartition    : Input does not satisfy the dedicated partition criteria
//                   and is handled by the general constraint machinery.
//
// ****************************************************************************

enum class PartitionType {
    RepStdAll             = 0,
    RepNoZero             = 1,
    RepShort              = 2,
    RepCapped             = 3,
    DistinctStdAll        = 4,
    DistinctMZ            = 5,
    DistinctOneZero       = 6,
    DistinctNoZero        = 7,
    DistinctCapped        = 8,
    DistinctCappedMZ      = 9,
    LengthOne             = 10,
    Multiset              = 11,
    CoarseGrained         = 12,
    CompRepNoZero         = 13,
    CompRepWeak           = 14,
    CompRepZero           = 15,
    CompDistinctNoZero    = 16,
    CompDistinctZero      = 17,
    CompDistinctWeak      = 18,
    CompDistinctMZWeak    = 19,
    CompDistinctCapped    = 20,
    CompDistinctCapWeak   = 21,
    CompDistinctCapMZ     = 22,
    CompDistinctCapMZWeak = 23,
    CompMultiset          = 24,
    PrmRepPartNoZ         = 25,
    PrmRepPart            = 26,
    PrmRepCapped          = 27,
    PrmDstPartNoZ         = 28,
    PrmDstPrtOneZ         = 29,
    PrmDstPartMZ          = 30,
    PrmDstPrtCap          = 31,
    PrmDstPrtCapMZ        = 32,
    PrmMultiset           = 33,
    NotMapped             = 34,
    NoSolution            = 35,
    NotPartition          = 36,

    // Appended to keep existing numeric values stable
    CompRepCapped         = 37,
    CompRepCapZero        = 38,
    CompRepWeakCap        = 39,
    CompMultisetZero      = 40,
    CompMultisetWeak      = 41,

    NumTypes              = 42
};

constexpr std::array<
    const char*,
    static_cast<size_t>(PartitionType::NumTypes)
> PTypeNames {{
    "RepStdAll",
    "RepNoZero",
    "RepShort",
    "RepCapped",
    "DistinctStdAll",
    "DistinctMZ",
    "DistinctOneZero",
    "DistinctNoZero",
    "DistinctCapped",
    "DistinctCappedMZ",
    "LengthOne",
    "Multiset",
    "CoarseGrained",
    "CompRepNoZero",
    "CompRepWeak",
    "CompRepZero",
    "CompDistinctNoZero",
    "CompDistinctZero",
    "CompDistinctWeak",
    "CompDistinctMZWeak",
    "CompDistinctCapped",
    "CompDistinctCapWeak",
    "CompDistinctCapMZ",
    "CompDistinctCapMZWeak",
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
    "CompRepCapped",
    "CompRepCapZero",
    "CompRepWeakCap",
    "CompMultisetZero",
    "CompMultisetWeak"
}};

const std::array<PartitionType, 4> NoCountAlgoPTypeArr{{
    PartitionType::NotMapped, PartitionType::NotPartition,
    PartitionType::NoSolution, PartitionType::CoarseGrained
}};

const std::array<PartitionType, 9> NoRankAlgoPTypeArr{{
    PartitionType::NotMapped, PartitionType::NotPartition,
    PartitionType::NoSolution, PartitionType::CoarseGrained,
    PartitionType::CompMultiset, PartitionType::CompMultisetWeak,
    PartitionType::CompMultisetZero, PartitionType::PrmMultiset,
    PartitionType::Multiset
}};

const std::array<PartitionType, 13> CappedPTypeArr{{
    PartitionType::RepCapped, PartitionType::DistinctCapped,
    PartitionType::DistinctCappedMZ, PartitionType::PrmRepCapped,
    PartitionType::PrmDstPrtCap, PartitionType::PrmDstPrtCapMZ,
    PartitionType::CompDistinctCapped, PartitionType::CompDistinctCapWeak,
    PartitionType::CompDistinctCapMZWeak, PartitionType::CompRepCapped,
    PartitionType::CompDistinctCapMZ, PartitionType::CompRepCapZero,
    PartitionType::CompRepWeakCap
}};

const std::array<PartitionType, 8> CompDistinctPTypeArr{{
    PartitionType::CompDistinctWeak, PartitionType::CompDistinctCapWeak,
    PartitionType::CompDistinctMZWeak, PartitionType::CompDistinctNoZero,
    PartitionType::CompDistinctZero, PartitionType::CompDistinctCapped,
    PartitionType::CompDistinctCapMZWeak, PartitionType::CompDistinctCapMZ
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
