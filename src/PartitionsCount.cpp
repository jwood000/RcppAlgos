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
#include <numeric>
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
            return std::make_unique<RepLenRstrctd>();
        } case PartitionType::DstctStdAll: {
            return std::make_unique<DistinctAll>();
        } case PartitionType::DstctMultiZero: {
            return std::make_unique<DistinctMZ>();
        } case PartitionType::DstctOneZero: {
            return std::make_unique<DistinctLen>();
        } case PartitionType::DstctNoZero: {
            return std::make_unique<DistinctLen>();
        } case PartitionType::DstctCapped: {
            return std::make_unique<DistinctLenRstrctd>();
        } case PartitionType::DstctCappedMZ: {
            return std::make_unique<DistinctRstrctdMZ>();
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
        } case PartitionType::CmpDstctCapped: {
            return std::make_unique<PermDstnctRstrctd>();
        } case PartitionType::CmpDstCapMZWeak: {
            return std::make_unique<PermDstnctRstrctdMZ>();
        } case PartitionType::CmpDstCapMZNotWk: {
            return std::make_unique<CompsDstnctRstrctdMZ>();
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
            return std::make_unique<PermDstnctRstrctd>();
        } case PartitionType::PrmDstPrtCapMZ: {
            return std::make_unique<PermDstnctRstrctdMZ>();
        } default: {
            return nullptr;
        }
    }
}

void CountClass::InitializeMpz() {
    if (size && width) {
        p2d.resize(width, std::vector<mpz_class>(size));
    } else if (size) {
        p1.resize(size);
        p2.resize(size);
    }
}

void DistinctLen::GetCount(mpz_class &res, int n, int m,
                           const std::vector<int> &allowed,
                           int strtLen, bool bLiteral) {

    if (cmp(res, 0) == 0 || cmp(res, Significand53) > 0) {
        CountPartsDistinctLen(res, p1, p2, n, m);
    } else {
        const double dblRes = CountPartsDistinctLen(n, m);
        res = dblRes;
    }
}

void DistinctLenRstrctd::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {

    if (cmp(res, 0) == 0 || cmp(res, Significand53) > 0) {
        CountPartsDistLenRstrctd(res, p2d, n, m, allowed);
    } else {
        res = CountPartsDistLenRstrctd(n, m, allowed);
    }
}

void DistinctMZ::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {

    if ((cmp(res, 0) == 0 || cmp(res, Significand53) > 0) && bLiteral) {
        CountPartsDistinctMultiZero(res, p1, p2, n, m, allowed, strtLen);
    } else if (cmp(res, 0) == 0 || cmp(res, Significand53) > 0) {
        CountPartsDistinctLen(res, p1, p2, n, m);
    } else if (bLiteral) {
        res = CountPartsDistinctMultiZero(n, m, allowed, strtLen);
    } else {
        res = CountPartsDistinctLen(n, m);
    }
}

void DistinctRstrctdMZ::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {

    if ((cmp(res, 0) == 0 || cmp(res, Significand53) > 0) && bLiteral) {
        CountPartsDistinctRstrctdMZ(res, p2d, n, m, allowed, strtLen);
    } else if (cmp(res, 0) == 0 || cmp(res, Significand53) > 0) {
        CountPartsDistLenRstrctd(res, p2d, n, m, allowed);
    } else if (bLiteral) {
        res = CountPartsDistinctRstrctdMZ(n, m, allowed, strtLen);
    } else {
        res = CountPartsDistLenRstrctd(n, m, allowed);
    }
}

void PermDstnctRstrctdMZ::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {

    if ((cmp(res, 0) == 0 || cmp(res, Significand53) > 0) && bLiteral) {
        CountPartsPermDistinctRstrctdMZ(res, p2d, n, m, allowed, strtLen);
    } else if (cmp(res, 0) == 0 || cmp(res, Significand53) > 0) {
        CountCompDistLenRstrctd(res, p2d, n, m, allowed);
    } else if (bLiteral) {
        res = CountPartsPermDistinctRstrctdMZ(n, m, allowed, strtLen);
    } else {
        res = CountCompDistLenRstrctd(n, m, allowed);
    }
}

void RepLen::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {

    if (cmp(res, 0) == 0 || cmp(res, Significand53) > 0) {
        CountPartsRepLen(res, p1, p2, n, m);
    } else {
        const double dblRes = CountPartsRepLen(n, m);
        res = dblRes;
    }
}

void RepLenRstrctd::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {

    if (cmp(res, 0) == 0 || cmp(res, Significand53) > 0) {
        CountPartsRepLenRstrctd(res, p2d, n, m, allowed);
    } else {
        const double dblRes = CountPartsRepLenRstrctd(n, m, allowed);
        res = dblRes;
    }
}

void DistinctAll::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {

    if (cmp(res, 0) == 0 || cmp(res, Significand53) > 0) {
        CountPartsDistinct(res, n, m);
    } else {
        const double dblRes = CountPartsDistinct(n, m);
        res = dblRes;
    }
}

void RepAll::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {

    if (cmp(res, 0) == 0 || cmp(res, Significand53) > 0) {
        CountPartsRep(res, n, m);
    } else {
        const double dblRes = CountPartsRep(n, m);
        res = dblRes;
    }
}

void PermDstnctRstrctd::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {

    if (cmp(res, 0) == 0 || cmp(res, Significand53) > 0) {
        CountCompDistLenRstrctd(res, p2d, n, m, allowed);
    } else {
        const double dblRes = CountCompDistLenRstrctd(n, m, allowed);
        res = dblRes;
    }
}

void CompsRepLen::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {
    CountCompsRepLen(res, n, m);
}

void CompsRepZero::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {

    if (bLiteral) {
        CountCompsRepZNotWk(res, n, m);
    } else {
        CountCompsRepLen(res, n, m);
    }
}

void CompsDistinctLen::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {
    CountCompsDistinctLen(res, p1, p2, n, m);
}

void CompsDistLenMZ::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {
    CountCompsDistinctMultiZero(res, p1, p2, n, m, allowed, strtLen);
}

void CompsDistLenMZWeak::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {
    CountCompsDistinctMZWeak(res, p1, p2, n, m, allowed, strtLen);
}

void CompsDstnctRstrctdMZ::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {

    if ((cmp(res, 0) == 0 || cmp(res, Significand53) > 0) && bLiteral) {
        CountCompsDistinctRstrctdMZ(res, p2d, n, m, allowed, strtLen);
    } else if (cmp(res, 0) == 0 || cmp(res, Significand53) > 0) {
        CountCompDistLenRstrctd(res, p2d, n, m, allowed);
    } else if (bLiteral) {
        res = CountCompsDistinctRstrctdMZ(n, m, allowed, strtLen);
    } else {
        res = CountCompDistLenRstrctd(n, m, allowed);
    }
}

double DistinctAll::GetCount(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {
    return CountPartsDistinct(n, m);
}

double DistinctLen::GetCount(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {
    return CountPartsDistinctLen(n, m);
}

double DistinctLenRstrctd::GetCount(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {
    return CountPartsDistLenRstrctd(n, m, allowed);
}

double DistinctMZ::GetCount(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {
    return CountPartsDistinctMultiZero(n, m, allowed, strtLen);
}

double DistinctRstrctdMZ::GetCount(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {
    return CountPartsDistinctRstrctdMZ(n, m, allowed, strtLen);
}

double RepAll::GetCount(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {
    return CountPartsRep(n, m);
}

double RepLen::GetCount(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {
    return CountPartsRepLen(n, m);
}

double RepLenRstrctd::GetCount(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {
    return CountPartsRepLenRstrctd(n, m, allowed);
}

double PermDstnctRstrctd::GetCount(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {
    return CountCompDistLenRstrctd(n, m, allowed);
}

double PermDstnctRstrctdMZ::GetCount(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {
    return CountPartsPermDistinctRstrctdMZ(n, m, allowed, strtLen);
}

double CompsRepLen::GetCount(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {
    return CountCompsRepLen(n, m);
}

double CompsRepZero::GetCount(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {
    return CountCompsRepZNotWk(n, m);
}

double CompsDistinctLen::GetCount(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {
    return CountCompsDistinctLen(n, m);
}

double CompsDistLenMZ::GetCount(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {
    return CountCompsDistinctMultiZero(n, m, allowed, strtLen);
}

double CompsDistLenMZWeak::GetCount(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {
    return CountCompsDistinctMZWeak(n, m, allowed, strtLen);
}

double CompsDstnctRstrctdMZ::GetCount(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {
    return CountCompsDistinctRstrctdMZ(n, m, allowed, strtLen);
}

bool IsTrueMultiset(PartitionType ptype) {

    switch(ptype) {
        case PartitionType::Multiset: {
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

void CountClass::SetArrSize(PartitionType ptype, int n, int m) {

    // N.B. We currently don't have a PrmRepCapped count function
    width = 0;
    size  = 0;

    switch (ptype) {
        case PartitionType::RepNoZero:
        case PartitionType::RepShort: {
            const int limit = std::min(n - m, m);
            CheckMultIsInt(2, m);
            CheckMultIsInt(2, limit);
            n = (n < 2 * m) ? 2 * limit : n;
            size = n + 1;
            return;
        }

        case PartitionType::DstctMultiZero:
        case PartitionType::DstctOneZero:
        case PartitionType::DstctNoZero: {
            CheckMultIsInt(1, n + 1);
            size = n + 1;
            return;
        }

        case PartitionType::RepCapped:
        case PartitionType::DstctCapped:
        case PartitionType::DstctCappedMZ:
        case PartitionType::CmpDstctCapped:
        case PartitionType::CmpDstCapMZWeak:
        case PartitionType::CmpDstCapMZNotWk:
        case PartitionType::PrmDstPrtCap:
        case PartitionType::PrmDstPrtCapMZ: {
            size  = n + 1;
            width = m + 1;
            return;
        }

        case PartitionType::CmpDstctNoZero:
        case PartitionType::CmpDstctMZWeak:
        case PartitionType::CmpDstctZNotWk: {
            CheckMultIsInt(1, n + 1);
            size = n + 1;
            return;
        }

        default: {
            width = 0;
            size  = 0;
            return;
        }
    }
}

int PartitionsCount(const std::vector<int> &Reps,
                    PartDesign &part, int lenV, bool bIsCount) {

    part.count = 0.0;
    part.numUnknown = false;
    part.bigCount = 0;

    if (part.ptype == PartitionType::NoSolution) {
        return 1;
    }

    const int strtLen = std::count_if(
        part.startZ.cbegin(), part.startZ.cend(), [](int i){return i > 0;}
    );

    // Returns false for all true multiset cases. Sometimes, part.isMult = true,
    // but it is really a distinct case with leading zeros.
    const bool bWorthIt = IsTrueMultiset(part.ptype);

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
        part.count = part.solnExist ?
            CountPartsMultiset(Reps, part.startZ, true, true) : 0;
        return 1;
    }

    std::unique_ptr<CountClass> Counter = MakeCount(part.ptype);

    if (Counter) {
        std::vector<int> allowed(part.cap);
        std::iota(allowed.begin(), allowed.end(), 1);

        part.count = Counter->GetCount(part.mapTar, part.width,
                                       allowed, strtLen);

        if (part.count > Significand53) {
            part.isGmp = true;
            Counter->SetArrSize(part.ptype, part.mapTar, part.width);
            Counter->InitializeMpz();
            Counter->GetCount(
                part.bigCount, part.mapTar, part.width, allowed, strtLen
            );
        }

        return 1;
    }

    // This shouldn't happen
    return -1;
}
