#include "Partitions/PartitionsCountMultiset.h"
#include "Partitions/PartitionsCountDistinct.h"
#include "Partitions/PartitionsCountSection.h"
#include "Partitions/BigPartsCountDistinct.h"
#include "Partitions/PartitionsCountRep.h"
#include "Partitions/BigPartsCountRep.h"
#include "Partitions/PartitionsCount.h"
#include "Permutations/PermuteCount.h"
#include "Combinations/ComboCount.h"
#include "CppConvert/Constants.h"  // Significand53
#include <algorithm>               // std::count_if, std::find
#include <memory>

std::unique_ptr<CountClass> MakeCount(PartitionType ptype) {

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
        } case PartitionType::CompRepNoZero: {
            return std::make_unique<CompsRepLen>();
        } case PartitionType::CompRepWeak: {
            return std::make_unique<CompsRepLen>();
        } case PartitionType::CmpRpZroNotWk: {
            return std::make_unique<CompsRepZero>();
        } case PartitionType::CmpDstctNoZero: {
            return std::make_unique<CompsDistinctLen>();
        } case PartitionType::CmpDstctZNotWk: {
            return std::make_unique<CompsDistLenMZ>();
        } case PartitionType::CmpDstctMZWeak: {
            return std::make_unique<CompsDistLenMZWeak>();
        } case PartitionType::PrmRepPart: {
            return std::make_unique<CompsRepLen>();
        } case PartitionType::PrmRepPartNoZ: {
            return std::make_unique<CompsRepLen>();
        } case PartitionType::PrmDstPartNoZ: {
            return std::make_unique<CompsDistinctLen>();
        } case PartitionType::PrmDstPrtOneZ: {
            return std::make_unique<CompsDistinctLen>();
        } case PartitionType::PrmDstPartMZ: {
            return std::make_unique<CompsDistLenMZWeak>();
        } case PartitionType::PrmDstPrtCap: {
            return std::make_unique<PermDstnctCap>();
        } case PartitionType::PrmDstPrtCapMZ: {
            return std::make_unique<PermDstnctCapMZ>();
        } default: {
            return nullptr;
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
        CountPartsDistinctLen(res, p1, p2, n, m);
    } else {
        const double dblRes = CountPartsDistinctLen(n, m);
        res = dblRes;
    }
}

void DistinctLenCap::GetCount(mpz_class &res, int n, int m, int cap,
                              int strtLen, bool bLiteral) {

    if (cmp(res, 0) == 0 || cmp(res, Significand53) > 0) {
        CountPartsDistinctLenCap(res, p1, p2, n, m, cap);
    } else {
        res = CountPartsDistinctLenCap(n, m, cap);
    }
}

void DistinctMZ::GetCount(mpz_class &res, int n, int m, int cap,
                          int strtLen, bool bLiteral) {

    if ((cmp(res, 0) == 0 || cmp(res, Significand53) > 0) && bLiteral) {
        CountPartsDistinctMultiZero(res, p1, p2, n, m, cap, strtLen);
    } else if (cmp(res, 0) == 0 || cmp(res, Significand53) > 0) {
        CountPartsDistinctLen(res, p1, p2, n, m);
    } else if (bLiteral) {
        res = CountPartsDistinctMultiZero(n, m, cap, strtLen);
    } else {
        res = CountPartsDistinctLen(n, m);
    }
}

void DistinctCapMZ::GetCount(mpz_class &res, int n, int m, int cap,
                             int strtLen, bool bLiteral) {

    if ((cmp(res, 0) == 0 || cmp(res, Significand53) > 0) && bLiteral) {
        CountPartsDistinctCapMZ(res, p1, p2, n, m, cap, strtLen);
    } else if (cmp(res, 0) == 0 || cmp(res, Significand53) > 0) {
        CountPartsDistinctLenCap(res, p1, p2, n, m, cap);
    } else if (bLiteral) {
        res = CountPartsDistinctCapMZ(n, m, cap, strtLen);
    } else {
        res = CountPartsDistinctLenCap(n, m, cap);
    }
}

void PermDstnctCapMZ::GetCount(mpz_class &res, int n, int m, int cap,
                               int strtLen, bool bLiteral) {

    if ((cmp(res, 0) == 0 || cmp(res, Significand53) > 0) && bLiteral) {
        CountPartsPermDistinctCapMZ(res, p1, p2, n, m, cap, strtLen);
    } else if (cmp(res, 0) == 0 || cmp(res, Significand53) > 0) {
        CountPartsPermDistinctCap(res, p1, p2, n, m, cap);
    } else if (bLiteral) {
        res = CountPartsPermDistinctCapMZ(n, m, cap, strtLen);
    } else {
        res = CountPartsPermDistinctCap(n, m, cap);
    }
}

void RepLen::GetCount(mpz_class &res, int n, int m, int cap,
                      int strtLen, bool bLiteral) {

    if (cmp(res, 0) == 0 || cmp(res, Significand53) > 0) {
        CountPartsRepLen(res, p1, p2, n, m);
    } else {
        const double dblRes = CountPartsRepLen(n, m);
        res = dblRes;
    }
}

void RepLenCap::GetCount(mpz_class &res, int n, int m, int cap,
                         int strtLen, bool bLiteral) {

    if (cmp(res, 0) == 0 || cmp(res, Significand53) > 0) {
        CountPartsRepLenCap(res, p1, p2, n, m, cap);
    } else {
        const double dblRes = CountPartsRepLenCap(n, m, cap);
        res = dblRes;
    }
}

void DistinctAll::GetCount(mpz_class &res, int n, int m, int cap,
                           int strtLen, bool bLiteral) {

    if (cmp(res, 0) == 0 || cmp(res, Significand53) > 0) {
        CountPartsDistinct(res, n, m);
    } else {
        const double dblRes = CountPartsDistinct(n, m);
        res = dblRes;
    }
}

void RepAll::GetCount(mpz_class &res, int n, int m, int cap,
                      int strtLen, bool bLiteral) {

    if (cmp(res, 0) == 0 || cmp(res, Significand53) > 0) {
        CountPartsRep(res, n, m);
    } else {
        const double dblRes = CountPartsRep(n, m);
        res = dblRes;
    }
}

void PermDstnctCap::GetCount(mpz_class &res, int n, int m, int cap,
                             int strtLen, bool bLiteral) {

    if (cmp(res, 0) == 0 || cmp(res, Significand53) > 0) {
        CountPartsPermDistinctCap(res, p1, p2, n, m, cap);
    } else {
        const double dblRes = CountPartsPermDistinctCap(n, m, cap);
        res = dblRes;
    }
}

void CompsRepLen::GetCount(mpz_class &res, int n, int m, int cap,
                           int strtLen, bool bLiteral) {
    CountCompsRepLen(res, n, m);
}

void CompsRepZero::GetCount(mpz_class &res, int n, int m, int cap,
                            int strtLen, bool bLiteral) {
    if (bLiteral) {
        CountCompsRepZNotWk(res, n, m);
    } else {
        CountCompsRepLen(res, n, m);
    }
}

void CompsDistinctLen::GetCount(mpz_class &res, int n, int m, int cap,
                                int strtLen, bool bLiteral) {
    CountCompsDistinctLen(res, p1, p2, n, m);
}

void CompsDistLenMZ::GetCount(mpz_class &res, int n, int m, int cap,
                                int strtLen, bool bLiteral) {
    CountCompsDistinctMultiZero(res, p1, p2, n, m, cap, strtLen);
}

void CompsDistLenMZWeak::GetCount(mpz_class &res, int n, int m, int cap,
                              int strtLen, bool bLiteral) {
    CountCompsDistinctMZWeak(res, p1, p2, n, m, cap, strtLen);
}

double DistinctAll::GetCount(int n, int m, int cap, int strtLen) {
    return CountPartsDistinct(n, m);
}

double DistinctLen::GetCount(int n, int m, int cap, int strtLen) {
    return CountPartsDistinctLen(n, m);
}

double DistinctLenCap::GetCount(int n, int m, int cap, int strtLen) {
    return CountPartsDistinctLenCap(n, m, cap);
}

double DistinctMZ::GetCount(int n, int m, int cap, int strtLen) {
    return CountPartsDistinctMultiZero(n, m, cap, strtLen);
}

double DistinctCapMZ::GetCount(int n, int m, int cap, int strtLen) {
    return CountPartsDistinctCapMZ(n, m, cap, strtLen);
}

double RepAll::GetCount(int n, int m, int cap, int strtLen) {
    return CountPartsRep(n, m);
}

double RepLen::GetCount(int n, int m, int cap, int strtLen) {
    return CountPartsRepLen(n, m);
}

double RepLenCap::GetCount(int n, int m, int cap, int strtLen) {
    return CountPartsRepLenCap(n, m, cap);
}

double PermDstnctCap::GetCount(int n, int m, int cap, int strtLen) {
    return CountPartsPermDistinctCap(n, m, cap);
}

double PermDstnctCapMZ::GetCount(int n, int m, int cap, int strtLen) {
    return CountPartsPermDistinctCapMZ(n, m, cap, strtLen);
}

double CompsRepLen::GetCount(int n, int m, int cap, int strtLen) {
    return CountCompsRepLen(n, m);
}

double CompsRepZero::GetCount(int n, int m, int cap, int strtLen) {
    return CountCompsRepZNotWk(n, m);
}

double CompsDistinctLen::GetCount(int n, int m, int cap, int strtLen) {
    return CountCompsDistinctLen(n, m);
}

double CompsDistLenMZ::GetCount(int n, int m, int cap, int strtLen) {
    return CountCompsDistinctMultiZero(n, m, cap, strtLen);
}

double CompsDistLenMZWeak::GetCount(int n, int m, int cap, int strtLen) {
    return CountCompsDistinctMZWeak(n, m, cap, strtLen);
}

bool OverTheBar(PartitionType ptype, double capNumIters, int n, int m) {

    // N.B. We currently don't have a PrmRepCapped count function
    constexpr double cutOff = 3.0;

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
        } case PartitionType::PrmDstPrtCap: {
            const double theBar = nChooseK(n, m);
            return (theBar / capNumIters) > cutOff;
        } case PartitionType::PrmDstPrtCapMZ: {
            const double theBar = nChooseK(n, m);
            return (theBar / capNumIters) > cutOff;
        } case PartitionType::Multiset: {
            return false;
        } case PartitionType::CompMultiset: {
            return false;
        } case PartitionType::PrmMultiset: {
            return false;
        } default: {
            return true;
        }
    }
}

void CountClass::SetArrSize(PartitionType ptype, int n, int m, int cap) {

    // N.B. We currently don't have a PrmRepCapped count function

    switch (ptype) {
        case PartitionType::RepNoZero: {
            const int limit = std::min(n - m, m);
            CheckMultIsInt(2, m);
            CheckMultIsInt(2, limit);
            n = (n < 2 * m) ? 2 * limit : n;
            size = n + 1;
            break;
        } case PartitionType::RepShort: {
            const int limit = std::min(n - m, m);
            CheckMultIsInt(2, m);
            CheckMultIsInt(2, limit);
            n = (n < 2 * m) ? 2 * limit : n;
            size = n + 1;
            break;
        } case PartitionType::RepCapped: {
            CheckMultIsInt(cap + 1, n + 1);
            size = (cap + 1) * (n + 1);
            break;
        } case PartitionType::DstctMultiZero: {
            CheckMultIsInt(1, n + 1);
            size = n + 1;
            break;
        } case PartitionType::DstctOneZero: {
            CheckMultIsInt(1, n + 1);
            size = n + 1;
            break;
        } case PartitionType::DstctNoZero: {
            CheckMultIsInt(1, n + 1);
            size = n + 1;
            break;
        } case PartitionType::DstctCapped: {
            CheckMultIsInt(cap + 1, n + 1);
            size = (cap + 1) * (n + 1);
            break;
        } case PartitionType::DstctCappedMZ: {
            CheckMultIsInt(cap + 1, n + 1);
            size = (cap + 1) * (n + 1);
            break;
        } case PartitionType::PrmDstPrtCap: {
            CheckMultIsInt(cap + 1, n + 1);
            size = (cap + 1) * (n + 1);
            break;
        } case PartitionType::PrmDstPrtCapMZ: {
            CheckMultIsInt(cap + 1, n + 1);
            size = (cap + 1) * (n + 1);
            break;
        } case PartitionType::CmpDstctNoZero: {
            CheckMultIsInt(1, n + 1);
            size = n + 1;
            break;
        } case PartitionType::CmpDstctMZWeak: {
            CheckMultIsInt(1, n + 1);
            size = n + 1;
            break;
        } case PartitionType::CmpDstctZNotWk: {
            CheckMultIsInt(1, n + 1);
            size = n + 1;
            break;
        } default: {
            size = 0;
            break;
        }
    }
}

int PartitionsCount(const std::vector<int> &Reps,
                    PartDesign &part, int lenV, bool bIsCount) {

    part.count = 0.0;
    part.numUnknown = false;
    part.bigCount = 0;

    const double capNumIters = static_cast<double>(part.mapTar + 1) *
                               static_cast<double>(part.width - 1) *
                               static_cast<double>(lenV + 1);

    const int strtLen = std::count_if(
        part.startZ.cbegin(), part.startZ.cend(), [](int i){return i > 0;}
    );

    // Returns false for all multiset cases and true for all non-capped cases
    const bool bWorthIt = OverTheBar(
        part.ptype, capNumIters, lenV, part.width
    );

    const auto no_algo_it = std::find(
        NoCountAlgoPTypeArr.cbegin(), NoCountAlgoPTypeArr.cend(), part.ptype
    );

    const bool forceCnt = (bIsCount || bWorthIt);

    if (!forceCnt || no_algo_it != NoCountAlgoPTypeArr.end()) {
        part.numUnknown = true;
        return 0;
    }

    if (part.ptype == PartitionType::LengthOne) {
        part.count = static_cast<int>(part.solnExist);
        return 1;
    } else if (part.ptype == PartitionType::Multiset) {
        part.count = part.solnExist ?
            CountPartsMultiset(Reps, part.startZ) : 0;
        return 1;
    } else if (part.ptype == PartitionType::CompMultiset) {
        // N.B. With PartitionType::Multiset we use the default
        // IsComp = false. Below, we set IsComp = true
        part.count = part.solnExist ?
            CountPartsMultiset(Reps, part.startZ, true, part.isWeak) : 0;
        return 1;
    } else if (part.ptype == PartitionType::PrmMultiset) {
        // See note above under CompMultiset
        part.count = (part.solnExist) ?
            CountPartsMultiset(Reps, part.startZ, true, true) : 0;
        return 1;
    }

    std::unique_ptr<CountClass> Counter = MakeCount(part.ptype);

    if (Counter) {
        part.count = Counter->GetCount(part.mapTar, part.width,
                                       part.cap, strtLen);

        if (part.count > Significand53) {
            part.isGmp = true;
            Counter->SetArrSize(
                part.ptype, part.mapTar, part.width, part.cap
            );
            Counter->InitializeMpz();
            Counter->GetCount(
                part.bigCount, part.mapTar, part.width, part.cap, strtLen
            );
        }

        return 1;
    }

    // This shouldn't happen
    return -1;
}
