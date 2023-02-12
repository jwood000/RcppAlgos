#include "Partitions/PartitionsCountMultiset.h"
#include "Partitions/PartitionsCountDistinct.h"
#include "Partitions/PartitionsCountSection.h"
#include "Partitions/BigPartsCountDistinct.h"
#include "Partitions/PartitionsCountRep.h"
#include "Partitions/BigPartsCountRep.h"
#include "Partitions/PartitionsCount.h"
#include "Permutations/PermuteCount.h"
#include "Combinations/ComboCount.h"
#include <memory>
#include "CppConvert/Constants.h"  // Significand53
#include <algorithm>               // std::count_if, std::find

constexpr double cutOff = 3.0;

std::unique_ptr<CountClass> MakeCount(PartitionType ptype, bool isComp) {

    if (isComp) {
        switch (ptype) {
            case PartitionType::RepStdAll: {
                return std::make_unique<CompsRepZero>();
            } case PartitionType::RepNoZero: {
                return std::make_unique<CompsRepLen>();
            } case PartitionType::RepShort: {
                return std::make_unique<CompsRepZero>();
            } default: {
                return nullptr;
            }
        }
    }

    switch (ptype) {
        case PartitionType::RepStdAll: {
            return std::make_unique<RepAll>();
        } case PartitionType::RepNoZero: {
            return std::make_unique<RepLen>();
        } case PartitionType::RepShort: {
            return std::make_unique<RepLen>();
        } case PartitionType::RepCapped: {
            return std::make_unique<RepLenCap>();
        } case PartitionType::DstctStdAll: {
            return std::make_unique<DistinctAll>();
        } case PartitionType::DstctMultiZero: {
            return std::make_unique<DistinctMZ>();
        } case PartitionType::DstctOneZero: {
            return std::make_unique<DistinctLen>();
        } case PartitionType::DstctNoZero: {
            return std::make_unique<DistinctLen>();
        } case PartitionType::DstctCapped: {
            return std::make_unique<DistinctLenCap>();
        } case PartitionType::DstctCappedMZ: {
            return std::make_unique<DistinctCapMZ>();
        } default: {
            return std::make_unique<RepAll>();
        }
    }
}


void CountClass::InitializeMpz() {
    if (size) {
        p1.resize(size);
        p2.resize(size);
    }
}

void DistinctLen::GetCount(mpz_class &res, int n, int m, int cap,
                           int strtLen, bool bLiteral) {

    if (cmp(res, 0) == 0 || cmp(res, Significand53) > 0) {
        CountPartsDistinctLen(res, p1, p2, n, m, cap, strtLen);
    } else {
        const double dblRes = CountPartsDistinctLen(n, m, cap, strtLen);
        res = dblRes;
    }
}

void DistinctLenCap::GetCount(mpz_class &res, int n, int m, int cap,
                              int strtLen, bool bLiteral) {

    if (cmp(res, 0) == 0 || cmp(res, Significand53) > 0) {
        CountPartsDistinctLenCap(res, p1, p2,
                                 n, m, cap, strtLen);
    } else {
        const double dblRes = CountPartsDistinctLenCap(n, m, cap, strtLen);
        res = dblRes;
    }
}

void DistinctMZ::GetCount(mpz_class &res, int n, int m, int cap,
                          int strtLen, bool bLiteral) {

    if (cmp(res, 0) == 0 || cmp(res, Significand53) > 0) {
        if (bLiteral) {
            CountPartsDistinctMultiZero(res, p1, p2,
                                        n, m, cap, strtLen);
        } else {
            CountPartsDistinctLen(res, p1, p2,
                                  n, m, cap, strtLen);
        }
    } else {
        const double dblRes = bLiteral ?
            CountPartsDistinctMultiZero(n, m, cap, strtLen) :
            CountPartsDistinctLen(n, m, cap, strtLen);
        res = dblRes;
    }
}

void DistinctCapMZ::GetCount(mpz_class &res, int n, int m, int cap,
                             int strtLen, bool bLiteral) {

    if (cmp(res, 0) == 0 || cmp(res, Significand53) > 0) {
        if (bLiteral) {
            CountPartsDistinctCapMZ(res, p1, p2,
                                    n, m, cap, strtLen);
        } else {
            CountPartsDistinctLenCap(res, p1, p2,
                                     n, m, cap, strtLen);
        }
    } else {
        const double dblRes = bLiteral ?
            CountPartsDistinctCapMZ(n, m, cap, strtLen) :
            CountPartsDistinctLenCap(n, m, cap, strtLen);
        res = dblRes;
    }
}

void RepLen::GetCount(mpz_class &res, int n, int m, int cap,
                      int strtLen, bool bLiteral) {

    if (cmp(res, 0) == 0 || cmp(res, Significand53) > 0) {
        CountPartsRepLen(res, p1, p2, n, m, cap, strtLen);
    } else {
        const double dblRes = CountPartsRepLen(n, m, cap, strtLen);
        res = dblRes;
    }
}

void RepLenCap::GetCount(mpz_class &res, int n, int m, int cap,
                         int strtLen, bool bLiteral) {

    if (cmp(res, 0) == 0 || cmp(res, Significand53) > 0) {
        CountPartsRepLenCap(res, p1, p2, n, m, cap, strtLen);
    } else {
        const double dblRes = CountPartsRepLenCap(n, m, cap, strtLen);
        res = dblRes;
    }
}

void DistinctAll::GetCount(mpz_class &res, int n, int m, int cap,
                           int strtLen, bool bLiteral) {

    if (cmp(res, 0) == 0 || cmp(res, Significand53) > 0) {
        CountPartsDistinct(res, n, m, cap, strtLen);
    } else {
        const double dblRes = CountPartsDistinct(n, m, cap, strtLen);
        res = dblRes;
    }
}

void RepAll::GetCount(mpz_class &res, int n, int m, int cap,
                      int strtLen, bool bLiteral) {

    if (cmp(res, 0) == 0 || cmp(res, Significand53) > 0) {
        CountPartsRep(res, n, m, cap, strtLen);
    } else {
        const double dblRes = CountPartsRep(n, m, cap, strtLen);
        res = dblRes;
    }
}

void CompsRepLen::GetCount(mpz_class &res, int n, int m, int cap,
                           int strtLen, bool bLiteral) {
    CountCompsRepLen(res, n, m, cap, strtLen);
}

void CompsRepZero::GetCount(mpz_class &res, int n, int m, int cap,
                            int strtLen, bool bLiteral) {
    if (bLiteral) {
        CountCompsRepZero(res, n, m, cap, strtLen);
    } else {
        CountCompsRepLen(res, n, m, cap, strtLen);
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

double CompsRepLen::GetCount(int n, int m, int cap, int strtLen) {
    return CountCompsRepLen(n, m, cap, strtLen);
}

double CompsRepZero::GetCount(int n, int m, int cap, int strtLen) {
    return CountCompsRepZero(n, m, cap, strtLen);
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
            CheckMultIsInt(2, m);
            CheckMultIsInt(2, limit);
            n = (n < 2 * m) ? 2 * limit : n;
            size = (n + 1);
            break;
        } case PartitionType::RepShort: {
            const int limit = std::min(n - m, m);
            CheckMultIsInt(2, m);
            CheckMultIsInt(2, limit);
            n = (n < 2 * m) ? 2 * limit : n;
            size = (n + 1);
            break;
        } case PartitionType::RepCapped: {
            CheckMultIsInt(cap + 1, n + 1);
            size = (cap + 1) * (n + 1);
            break;
        } case PartitionType::DstctMultiZero: {
            CheckMultIsInt(1, n + 1);
            size = (n + 1);
            break;
        } case PartitionType::DstctOneZero: {
            CheckMultIsInt(1, n + 1);
            size = (n + 1);
            break;
        } case PartitionType::DstctNoZero: {
            CheckMultIsInt(1, n + 1);
            size = (n + 1);
            break;
        } case PartitionType::DstctCapped: {
            CheckMultIsInt(cap + 1, n + 1);
            size = (cap + 1) * (n + 1);
            break;
        } case PartitionType::DstctCappedMZ: {
            CheckMultIsInt(cap + 1, n + 1);
            size = (cap + 1) * (n + 1);
            break;
        } default: {
            size = 0;
            break;
        }
    }
}

void PartitionsCount(const std::vector<int> &Reps,
                     PartDesign &part, int lenV, bool bIsCount) {

    part.count = 0.0;
    part.numUnknown = false;
    part.bigCount = 0;

    const double capNumIters = static_cast<double>(part.mapTar + 1) *
                               static_cast<double>(part.width - 1) *
                               static_cast<double>(lenV + 1);

    const int strtLen = std::count_if(part.startZ.cbegin(),
                                      part.startZ.cend(),
                                      [](int i){return i > 0;});

    const bool bWorthIt = OverTheBar(part.ptype, capNumIters,
                                     lenV, part.width);

    if (part.ptype == PartitionType::LengthOne) {
        part.count = static_cast<int>(part.solnExist);
    } else if (part.isComb && part.ptype != PartitionType::Multiset) {
        if (bIsCount || bWorthIt) {
            std::unique_ptr<CountClass> myClass = MakeCount(part.ptype);
            part.count = myClass->GetCount(part.mapTar, part.width,
                                           part.cap, strtLen);

            if (part.count > Significand53) {
                part.isGmp = true;

                if (part.ptype != PartitionType::RepStdAll &&
                    part.ptype != PartitionType::DstctStdAll) {

                    myClass->SetArrSize(part.ptype, part.mapTar,
                                        part.width, part.cap);
                    myClass->InitializeMpz();
                    myClass->GetCount(part.bigCount, part.mapTar,
                                      part.width, part.cap, strtLen);
                } else {
                    myClass->GetCount(part.bigCount, part.mapTar,
                                      part.width, part.cap, strtLen);
                }
            }
        } else {
            part.numUnknown = true;
        }
    } else if (part.isComb) {
        if (bIsCount) {
            part.count = (part.solnExist) ?
                         CountPartsMultiset(Reps, part.startZ) : 0;
        } else {
            part.numUnknown = true;
        }
    } else {
        const auto it = std::find(DistPTypeArr.cbegin(),
                                  DistPTypeArr.cend(), part.ptype);

        if (part.isComp) {
            std::unique_ptr<CountClass> myClass = MakeCount(part.ptype, true);

            if (myClass) {
                part.count = myClass->GetCount(part.mapTar, part.width,
                                               part.cap, strtLen);

                if (part.count > Significand53) {
                    part.isGmp = true;
                    myClass->GetCount(part.bigCount, part.mapTar,
                                      part.width, part.cap, strtLen);
                }
            } else {
                part.numUnknown = true;
            }
        } else if (part.isRep && part.ptype != PartitionType::RepCapped) {
            part.count = CountCompsRepLen(
                part.mapTar + static_cast<int>(part.mapIncZero) * part.width,
                part.width, part.cap, strtLen
            );
        } else if (part.ptype == PartitionType::DstctCappedMZ ||
                   part.ptype == PartitionType::DstctCapped) {

            if (bIsCount || bWorthIt) {
                part.count = CountPartsPermDistinctCap(part.startZ, part.cap,
                                                       part.mapTar, part.width,
                                                       part.mapIncZero);
            } else {
                part.numUnknown = true;
            }
        } else if (it != DistPTypeArr.cend()) {
            part.count = CountPartsPermDistinct(part.startZ,
                                                part.mapTar, part.width,
                                                part.mapIncZero);
        } else {
            part.numUnknown = true;
        }
    }
}
