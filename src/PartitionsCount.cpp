#include "Partitions/PartitionsCountMultiset.h"
#include "Partitions/PartitionsCountDistinct.h"
#include "Partitions/BigPartsCountDistinct.h"
#include "Partitions/PartitionsCountRep.h"
#include "Partitions/BigPartsCountRep.h"
#include "Partitions/PartitionsCount.h"
#include "Permutations/PermuteCount.h"
#include "Combinations/ComboCount.h"
#include "CleanConvert.h"  // Significand53
#include <algorithm>       // std::count_if, std::find

constexpr double cutOff = 3.0;

std::unique_ptr<CountClass> CountClass::MakeCount(PartitionType ptype) {

    switch (ptype) {
        case PartitionType::RepStdAll: {
            return FromCpp14::make_unique<RepAll>();
        } case PartitionType::RepNoZero: {
            return FromCpp14::make_unique<RepLen>();
        } case PartitionType::RepShort: {
            return FromCpp14::make_unique<RepLen>();
        } case PartitionType::RepCapped: {
            return FromCpp14::make_unique<RepLenCap>();
        } case PartitionType::DstctStdAll: {
            return FromCpp14::make_unique<DistinctAll>();
        } case PartitionType::DstctMultiZero: {
            return FromCpp14::make_unique<DistinctMZ>();
        } case PartitionType::DstctOneZero: {
            return FromCpp14::make_unique<DistinctLen>();
        } case PartitionType::DstctNoZero: {
            return FromCpp14::make_unique<DistinctLen>();
        } case PartitionType::DstctCapped: {
            return FromCpp14::make_unique<DistinctLenCap>();
        } case PartitionType::DstctCappedMZ: {
            return FromCpp14::make_unique<DistinctCapMZ>();
        } default: {
            return FromCpp14::make_unique<RepAll>();
        }
    }
}

void CountClass::ClearMpz() {
    for (int i = 0; i < size; ++i) {
        mpz_clear(p1[i]);
        mpz_clear(p2[i]);
    }
}

void CountClass::InitializeMpz() {
    if (size) {
        p1 = FromCpp14::make_unique<mpz_t[]>(size);
        p2 = FromCpp14::make_unique<mpz_t[]>(size);

        for (int i = 0; i < size; ++i) {
            mpz_init(p1[i]);
            mpz_init(p2[i]);
        }
    }
}

void DistinctLen::GetCount(mpz_t res, int n, int m, int cap,
                           int strtLen, bool bLiteral) {

    if (mpz_cmp_ui(res, 0u) == 0 || mpz_cmp_d(res, Significand53) > 0) {
        CountPartsDistinctLen(res, p1.get(), p2.get(), n, m, cap, strtLen);
    } else {
        const double dblRes = CountPartsDistinctLen(n, m, cap, strtLen);
        mpz_set_d(res, dblRes);
    }
}

void DistinctLenCap::GetCount(mpz_t res, int n, int m, int cap,
                              int strtLen, bool bLiteral) {

    if (mpz_cmp_ui(res, 0u) == 0 || mpz_cmp_d(res, Significand53) > 0) {
        CountPartsDistinctLenCap(res, p1.get(), p2.get(),
                                 n, m, cap, strtLen);
    } else {
        const double dblRes = CountPartsDistinctLenCap(n, m, cap, strtLen);
        mpz_set_d(res, dblRes);
    }
}

void DistinctMZ::GetCount(mpz_t res, int n, int m, int cap,
                          int strtLen, bool bLiteral) {

    if (mpz_cmp_ui(res, 0u) == 0 || mpz_cmp_d(res, Significand53) > 0) {
        if (bLiteral) {
            CountPartsDistinctMultiZero(res, p1.get(), p2.get(),
                                        n, m, cap, strtLen);
        } else {
            CountPartsDistinctLen(res, p1.get(), p2.get(),
                                  n, m, cap, strtLen);
        }
    } else {
        const double dblRes = bLiteral ?
            CountPartsDistinctMultiZero(n, m, cap, strtLen) :
            CountPartsDistinctLen(n, m, cap, strtLen);
        mpz_set_d(res, dblRes);
    }
}

void DistinctCapMZ::GetCount(mpz_t res, int n, int m, int cap,
                             int strtLen, bool bLiteral) {

    if (mpz_cmp_ui(res, 0u) == 0 || mpz_cmp_d(res, Significand53) > 0) {
        if (bLiteral) {
            CountPartsDistinctCapMZ(res, p1.get(), p2.get(),
                                    n, m, cap, strtLen);
        } else {
            CountPartsDistinctLenCap(res, p1.get(), p2.get(),
                                     n, m, cap, strtLen);
        }
    } else {
        const double dblRes = bLiteral ?
            CountPartsDistinctCapMZ(n, m, cap, strtLen) :
            CountPartsDistinctLenCap(n, m, cap, strtLen);
        mpz_set_d(res, dblRes);
    }
}

void RepLen::GetCount(mpz_t res, int n, int m, int cap,
                      int strtLen, bool bLiteral) {

    if (mpz_cmp_ui(res, 0u) == 0 || mpz_cmp_d(res, Significand53) > 0) {
        CountPartsRepLen(res, p1.get(), p2.get(), n, m, cap, strtLen);
    } else {
        const double dblRes = CountPartsRepLen(n, m, cap, strtLen);
        mpz_set_d(res, dblRes);
    }
}

void RepLenCap::GetCount(mpz_t res, int n, int m, int cap,
                         int strtLen, bool bLiteral) {

    if (mpz_cmp_ui(res, 0u) == 0 || mpz_cmp_d(res, Significand53) > 0) {
        CountPartsRepLenCap(res, p1.get(), p2.get(), n, m, cap, strtLen);
    } else {
        const double dblRes = CountPartsRepLenCap(n, m, cap, strtLen);
        mpz_set_d(res, dblRes);
    }
}

void DistinctAll::GetCount(mpz_t res, int n, int m, int cap,
                           int strtLen, bool bLiteral) {

    if (mpz_cmp_ui(res, 0u) == 0 || mpz_cmp_d(res, Significand53) > 0) {
        CountPartsDistinct(res, n, m, cap, strtLen);
    } else {
        const double dblRes = CountPartsDistinct(n, m, cap, strtLen);
        mpz_set_d(res, dblRes);
    }
}

void RepAll::GetCount(mpz_t res, int n, int m, int cap,
                      int strtLen, bool bLiteral) {

    if (mpz_cmp_ui(res, 0u) == 0 || mpz_cmp_d(res, Significand53) > 0) {
        CountPartsRep(res, n, m, cap, strtLen);
    } else {
        const double dblRes = CountPartsRep(n, m, cap, strtLen);
        mpz_set_d(res, dblRes);
    }
}

double DistinctAll::GetCount(int n, int m, int cap, int strtLen) {
    return CountPartsDistinct(n, m, cap, strtLen);
}

double DistinctLen::GetCount(int n, int m, int cap, int strtLen) {
    return CountPartsDistinctLen(n, m, cap, strtLen);
}

double DistinctLenCap::GetCount(int n, int m, int cap, int strtLen) {
    return CountPartsDistinctLenCap(n, m, cap, strtLen);
}

double DistinctMZ::GetCount(int n, int m, int cap, int strtLen) {
    return CountPartsDistinctMultiZero(n, m, cap, strtLen);
}

double DistinctCapMZ::GetCount(int n, int m, int cap, int strtLen) {
    return CountPartsDistinctCapMZ(n, m, cap, strtLen);
}

double RepAll::GetCount(int n, int m, int cap, int strtLen) {
    return CountPartsRep(n, m, cap, strtLen);
}

double RepLen::GetCount(int n, int m, int cap, int strtLen) {
    return CountPartsRepLen(n, m, cap, strtLen);
}

double RepLenCap::GetCount(int n, int m, int cap, int strtLen) {
    return CountPartsRepLenCap(n, m, cap, strtLen);
}

bool OverTheBar(PartitionType ptype, double capNumIters, int n, int m) {
    switch(ptype) {
        case PartitionType::RepCapped: {
            const double theBar = NumCombsWithRep(n, m);
            return (theBar / capNumIters) > cutOff;
        } case PartitionType::DstctCapped: {
            const double theBar = nChooseK(n, m);
            return (theBar / capNumIters) > cutOff;
        } case PartitionType::DstctCappedMZ: {
            const double theBar = nChooseK(n, m);
            return (theBar / capNumIters) > cutOff;
        } default: {
            return true;
        }
    }
}

void CountClass::SetArrSize(PartitionType ptype, int n, int m, int cap) {

    switch (ptype) {
        case PartitionType::RepNoZero: {
            const int limit = std::min(n - m, m);
            n = (n < 2 * m) ? 2 * limit : n;
            size = (n + 1);
            break;
        } case PartitionType::RepShort: {
            const int limit = std::min(n - m, m);
            n = (n < 2 * m) ? 2 * limit : n;
            size = (n + 1);
            break;
        } case PartitionType::RepCapped: {
            size = (cap + 1) * (n + 1);
            break;
        } case PartitionType::DstctMultiZero: {
            size = (n + 1);
            break;
        } case PartitionType::DstctOneZero: {
            size = (n + 1);
            break;
        } case PartitionType::DstctNoZero: {
            size = (n + 1);
            break;
        } case PartitionType::DstctCapped: {
            size = (cap + 1) * (n + 1);
            break;
        } case PartitionType::DstctCappedMZ: {
            size = (cap + 1) * (n + 1);
            break;
        } default: {
            size = 0;
            break;
        }
    }
}

void PartitionsCount(const std::vector<int> &Reps, PartDesign &part,
                     int lenV, bool bCalcDifficult, bool IsComb) {

    part.count = 0.0;
    mpz_init(part.bigCount);
    const double capNumIters = static_cast<double>(part.mapTar + 1) *
                               static_cast<double>(part.width - 1) *
                               static_cast<double>(lenV + 1);

    const bool bWorthIt = OverTheBar(part.ptype, capNumIters, lenV, part.width);

    if (IsComb && part.ptype != PartitionType::Multiset) {
        if (bCalcDifficult || bWorthIt) {
            const int strtLen = std::count_if(part.startZ.cbegin(),
                                              part.startZ.cend(),
                                              [](int i){return i > 0;});

            const int cap = lenV - part.mapIncZero;

            CountClass init;
            std::unique_ptr<CountClass> myClass = init.MakeCount(part.ptype);
            part.count = myClass->GetCount(part.mapTar, part.width,
                                           cap, strtLen);

            if (part.count > Significand53) {
                part.isGmp = true;

                if (part.ptype != PartitionType::RepStdAll &&
                    part.ptype != PartitionType::DstctStdAll) {

                    myClass->SetArrSize(part.ptype, part.mapTar,
                                        part.width, cap);
                    myClass->InitializeMpz();
                    myClass->GetCount(part.bigCount, part.mapTar,
                                      part.width, cap, strtLen);
                    myClass->ClearMpz();
                } else {
                    myClass->GetCount(part.bigCount, part.mapTar,
                                      part.width, cap, strtLen);
                }
            }
        } else {
            part.numUnknown = true;
        }
    } else if (IsComb && part.ptype == PartitionType::Multiset) {
        if (bCalcDifficult && part.solnExist) {
            part.count = CountPartsMultiset(Reps, part.startZ);
        } else {
            part.numUnknown = true;
        }
    } else {
        const auto it = std::find(DistPTypeArr.cbegin(),
                                  DistPTypeArr.cend(), part.ptype);
        if (part.isRep) {
            if (part.ptype != PartitionType::RepCapped) {
                part.count = CountPartsPermRep(part.mapTar, part.width,
                                               part.mapIncZero);
            } else {
                part.numUnknown = true;
            }
        } else if (part.ptype == PartitionType::DstctCapped ||
            part.ptype == PartitionType::DstctCappedMZ) {

            if (bCalcDifficult || bWorthIt) {
                part.count = CountPartsPermDistinctCap(part.startZ,
                                                       lenV - part.mapIncZero,
                                                       part.mapTar, part.width,
                                                       part.mapIncZero);
            } else {
                part.numUnknown = true;
            }
        } else if (it != DistPTypeArr.cend()) {
            part.count = CountPartsPermDistinct(part.startZ, part.mapTar,
                                                part.width,
                                                part.mapIncZero);
        } else {
            part.numUnknown = true;
        }
    }
}
