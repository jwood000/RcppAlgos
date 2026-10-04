#include "Partitions/PartitionsCountMultiset.h"
#include "Partitions/PartitionsCountDistinct.h"
#include "Partitions/PartitionsCountSection.h"
#include "Partitions/BigPartsCountDistinct.h"
#include "Partitions/BigPartsCountMultiset.h"
#include "Partitions/PartitionsCountRep.h"
#include "Partitions/MultisetCountClass.h"
#include "Partitions/BigPartsCountRep.h"
#include "Partitions/PartitionsCount.h"
#include "Permutations/PermuteCount.h"
#include "Combinations/ComboCount.h"
#include "SetUpUtils.h"  // IsBeyondBound
#include <algorithm>     // std::count_if, std::find, std::min
#include <numeric>       // std::iota
#include <memory>        // std::make_unique, std::unique_ptr

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
        } case PartitionType::DistinctStdAll: {
            return std::make_unique<DistinctAll>();
        } case PartitionType::DistinctMZ: {
            return std::make_unique<DistinctMZ>();
        } case PartitionType::DistinctOneZero: {
            return std::make_unique<DistinctLen>();
        } case PartitionType::DistinctNoZero: {
            return std::make_unique<DistinctLen>();
        } case PartitionType::DistinctCapped: {
            return std::make_unique<DistinctLenRstrctd>();
        } case PartitionType::DistinctCappedMZ: {
            return std::make_unique<DistinctRstrctdMZ>();
        } case PartitionType::CompRepNoZero: {
            return std::make_unique<CompsRepLen>();
        } case PartitionType::CompRepWeak: {
            return std::make_unique<CompsRepLen>();
        } case PartitionType::CompRepCapped: {
            return std::make_unique<CompsRepLenCap>();
        } case PartitionType::CompRepWeakCap: {
            return std::make_unique<CompsRepLenCap>();
        } case PartitionType::CompRepCapZero: {
            return std::make_unique<CompsRepZeroCap>();
        } case PartitionType::CompRepZero: {
            return std::make_unique<CompsRepZero>();
        } case PartitionType::CompDistinctNoZero: {
            return std::make_unique<CompsDistinctLen>();
        } case PartitionType::CompDistinctZero: {
            return std::make_unique<CompsDistLenMZ>();
        } case PartitionType::CompDistinctWeak: {
            return std::make_unique<CompsDistinctLen>();
        } case PartitionType::CompDistinctMZWeak: {
            return std::make_unique<CompsDistLenMZWeak>();
        } case PartitionType::CompDistinctCapped: {
            return std::make_unique<PermDstnctRstrctd>();
        } case PartitionType::CompDistinctCapWeak: {
            return std::make_unique<PermDstnctRstrctd>();
        } case PartitionType::CompDistinctCapMZWeak: {
            return std::make_unique<PermDstnctRstrctdMZ>();
        } case PartitionType::CompDistinctCapMZ: {
            return std::make_unique<CompsDstnctRstrctdMZ>();
        } case PartitionType::PrmRepPart: {
            return std::make_unique<CompsRepLen>();
        } case PartitionType::PrmRepPartNoZ: {
            return std::make_unique<CompsRepLen>();
        } case PartitionType::PrmRepCapped: {
            return std::make_unique<CompsRepLenCap>();
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

std::unique_ptr<MultisetCountClass> MakeMultisetCount(
    PartitionType ptype, const std::vector<int>& reps
) {

    switch (ptype) {
        case PartitionType::Multiset: {
            return std::make_unique<PartsMultiset>(reps);
        } case PartitionType::CompMultiset: {
            return std::make_unique<CompsMultiset>(reps);
        } case PartitionType::CompMultisetZero: {
            return std::make_unique<CompsMultisetZero>(reps);
        } case PartitionType::CompMultisetWeak: {
            return std::make_unique<CompsMultisetWeak>(reps);
        } case PartitionType::PrmMultiset: {
            // PrmMultiset always treats zero as an actual part, so its count
            // semantics are equivalent to standard multiset compositions.
            return std::make_unique<CompsMultiset>(reps);
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

    double dblRes = 0;
    bool computedDouble = false;

    if (cmp(res, Significand53) < 0) {
        dblRes = CountPartsDistinctLen(n, m);
        computedDouble = true;
    }

    if (!computedDouble || IsBeyondBound(dblRes)) {
        CountPartsDistinctLen(res, p1, p2, n, m);
    } else {
        res = dblRes;
    }
}

void DistinctLenRstrctd::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {

    double dblRes = 0;
    bool computedDouble = false;

    if (cmp(res, Significand53) < 0) {
        dblRes = CountPartsDistLenRstrctd(n, m, allowed);
        computedDouble = true;
    }

    if (!computedDouble || IsBeyondBound(dblRes)) {
        CountPartsDistLenRstrctd(res, p2d, n, m, allowed);
    } else {
        res = dblRes;
    }
}

void DistinctMZ::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {

    double dblRes = 0;
    bool computedDouble = false;

    if (cmp(res, Significand53) < 0) {
        if (bLiteral) {
            dblRes = CountPartsDistinctMZ(n, m, allowed, strtLen);
        } else {
            dblRes = CountPartsDistinctLen(n, m);
        }

        computedDouble = true;
    }

    if (!computedDouble || IsBeyondBound(dblRes)) {
        if (bLiteral) {
            CountPartsDistinctMZ(res, p1, p2, n, m, allowed, strtLen);
        } else {
            CountPartsDistinctLen(res, p1, p2, n, m);
        }
    } else {
        res = dblRes;
    }
}

void DistinctRstrctdMZ::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {

    double dblRes = 0;
    bool computedDouble = false;

    if (cmp(res, Significand53) < 0) {
        if (bLiteral) {
            dblRes = CountPartsDistinctRstrctdMZ(n, m, allowed, strtLen);
        } else {
            dblRes = CountPartsDistLenRstrctd(n, m, allowed);
        }

        computedDouble = true;
    }

    if (!computedDouble || IsBeyondBound(dblRes)) {
        if (bLiteral) {
            CountPartsDistinctRstrctdMZ(res, p2d, n, m, allowed, strtLen);
        } else {
            CountPartsDistLenRstrctd(res, p2d, n, m, allowed);
        }
    } else {
        res = dblRes;
    }
}

void PermDstnctRstrctdMZ::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {

    double dblRes = 0;
    bool computedDouble = false;

    if (cmp(res, Significand53) < 0) {
        if (bLiteral) {
            dblRes = CountPartsPermDistinctRstrctdMZ(n, m, allowed, strtLen);
        } else {
            dblRes = CountCompsDistLenRstrctd(n, m, allowed);
        }

        computedDouble = true;
    }

    if (!computedDouble || IsBeyondBound(dblRes)) {
        if (bLiteral) {
            CountPartsPermDistinctRstrctdMZ(res, p2d, n, m, allowed, strtLen);
        } else {
            CountCompsDistLenRstrctd(res, p2d, n, m, allowed);
        }
    } else {
        res = dblRes;
    }
}

void RepLen::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {

    double dblRes = 0;
    bool computedDouble = false;

    if (cmp(res, Significand53) < 0) {
        dblRes = CountPartsRepLen(n, m);
        computedDouble = true;
    }

    if (!computedDouble || IsBeyondBound(dblRes)) {
        CountPartsRepLen(res, p1, p2, n, m);
    } else {
        res = dblRes;
    }
}

void RepLenRstrctd::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {

    double dblRes = 0;
    bool computedDouble = false;

    if (cmp(res, Significand53) < 0) {
        dblRes = CountPartsRepLenRstrctd(n, m, allowed);
        computedDouble = true;
    }

    if (!computedDouble || IsBeyondBound(dblRes)) {
        CountPartsRepLenRstrctd(res, p2d, n, m, allowed);
    } else {
        res = dblRes;
    }
}

void DistinctAll::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {

    double dblRes = 0;
    bool computedDouble = false;

    if (cmp(res, Significand53) < 0) {
        dblRes = CountPartsDistinct(n, m);
        computedDouble = true;
    }

    if (!computedDouble || IsBeyondBound(dblRes)) {
        CountPartsDistinct(res, n, m);
    } else {
        res = dblRes;
    }
}

void RepAll::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {

    double dblRes = 0;
    bool computedDouble = false;

    if (cmp(res, Significand53) < 0) {
        dblRes = CountPartsRep(n, m);
        computedDouble = true;
    }

    if (!computedDouble || IsBeyondBound(dblRes)) {
        CountPartsRep(res, n, m);
    } else {
        res = dblRes;
    }
}

void PartsMultiset::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {

    double dblRes = 0;
    bool computedDouble = false;

    if (cmp(res, Significand53) < 0) {
        dblRes = CountPartsMultiset(n, m, allowed, Reps);
        computedDouble = true;
    }

    if (!computedDouble || IsBeyondBound(dblRes)) {
        CountPartsMultiset(res, n, m, allowed, Reps);
    } else {
        res = dblRes;
    }
}

void PermDstnctRstrctd::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {

    double dblRes = 0;
    bool computedDouble = false;

    if (cmp(res, Significand53) < 0) {
        dblRes = CountCompsDistLenRstrctd(n, m, allowed);
        computedDouble = true;
    }

    if (!computedDouble || IsBeyondBound(dblRes)) {
        CountCompsDistLenRstrctd(res, p2d, n, m, allowed);
    } else {
        res = dblRes;
    }
}

void CompsRepLen::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {
    CountCompsRepLen(res, n, m);
}

void CompsRepLenCap::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {
    CountCompsRepLenCap(res, n, m, allowed);
}

void CompsRepZero::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {

    if (bLiteral) {
        CountCompsRepZero(res, n, m);
    } else {
        CountCompsRepLen(res, n, m);
    }
}

void CompsRepZeroCap::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {

    if (bLiteral) {
        CountCompsRepCapZero(res, n, m, allowed);
    } else {
        CountCompsRepLenCap(res, n, m, allowed);
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
    CountCompsDistinctMZ(res, p1, p2, n, m, allowed, strtLen);
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

    double dblRes = 0;
    bool computedDouble = false;

    if (cmp(res, Significand53) < 0) {
        if (bLiteral) {
            dblRes = CountCompsDistinctRstrctdMZ(n, m, allowed, strtLen);
        } else {
            dblRes = CountCompsDistLenRstrctd(n, m, allowed);
        }
        computedDouble = true;
    }

    if (!computedDouble || IsBeyondBound(dblRes)) {
        if (bLiteral) {
            CountCompsDistinctRstrctdMZ(res, p2d, n, m, allowed, strtLen);
        } else {
            CountCompsDistLenRstrctd(res, p2d, n, m, allowed);
        }
    } else {
        res = dblRes;
    }
}

void CompsMultiset::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {

    double dblRes = 0;
    bool computedDouble = false;

    if (cmp(res, Significand53) < 0) {
        dblRes = CountCompsMultiset(n, m, allowed, Reps);
        computedDouble = true;
    }

    if (!computedDouble || IsBeyondBound(dblRes)) {
        CountCompsMultiset(res, n, m, allowed, Reps);
    } else {
        res = dblRes;
    }
}

void CompsMultisetZero::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {

    double dblRes = 0;
    bool computedDouble = false;

    if (cmp(res, Significand53) < 0) {
        dblRes = CountCompsMultisetZero(n, m, allowed, Reps, strtLen);
        computedDouble = true;
    }

    if (!computedDouble || IsBeyondBound(dblRes)) {
        CountCompsMultisetZero(res, n, m, allowed, Reps, strtLen);
    } else {
        res = dblRes;
    }
}

void CompsMultisetWeak::GetCount(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen, bool bLiteral
) {

    double dblRes = 0;
    bool computedDouble = false;

    if (cmp(res, Significand53) < 0) {
        dblRes = CountCompsMultisetWeak(n, m, allowed, Reps, strtLen);
        computedDouble = true;
    }

    if (!computedDouble || IsBeyondBound(dblRes)) {
        CountCompsMultisetWeak(res, n, m, allowed, Reps, strtLen);
    } else {
        res = dblRes;
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
    return CountPartsDistinctMZ(n, m, allowed, strtLen);
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

double PartsMultiset::GetCount(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {
    return CountPartsMultiset(n, m, allowed, Reps);
}

double PermDstnctRstrctd::GetCount(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {
    return CountCompsDistLenRstrctd(n, m, allowed);
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

double CompsRepLenCap::GetCount(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {
    return CountCompsRepLenCap(n, m, allowed);
}

double CompsRepZero::GetCount(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {
    return CountCompsRepZero(n, m);
}

double CompsRepZeroCap::GetCount(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {
    return CountCompsRepCapZero(n, m, allowed);
}

double CompsDistinctLen::GetCount(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {
    return CountCompsDistinctLen(n, m);
}

double CompsDistLenMZ::GetCount(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {
    return CountCompsDistinctMZ(n, m, allowed, strtLen);
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

double CompsMultiset::GetCount(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {
    return CountCompsMultiset(n, m, allowed, Reps);
}

double CompsMultisetZero::GetCount(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {
    return CountCompsMultisetZero(n, m, allowed, Reps, strtLen);
}

double CompsMultisetWeak::GetCount(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {
    return CountCompsMultisetWeak(n, m, allowed, Reps, strtLen);
}

void CountClass::SetArrSize(PartitionType ptype, int n, int m) {

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

        case PartitionType::DistinctMZ:
        case PartitionType::DistinctOneZero:
        case PartitionType::DistinctNoZero: {
            CheckMultIsInt(1, n + 1);
            size = n + 1;
            return;
        }

        case PartitionType::RepCapped:
        case PartitionType::DistinctCapped:
        case PartitionType::DistinctCappedMZ:
        case PartitionType::CompDistinctCapped:
        case PartitionType::CompDistinctCapWeak:
        case PartitionType::CompRepCapped:
        case PartitionType::CompDistinctCapMZWeak:
        case PartitionType::CompDistinctCapMZ:
        case PartitionType::PrmDstPrtCap:
        case PartitionType::PrmRepCapped:
        case PartitionType::PrmDstPrtCapMZ: {
            size  = n + 1;
            width = m + 1;
            return;
        }

        case PartitionType::CompDistinctWeak:
        case PartitionType::CompDistinctNoZero:
        case PartitionType::CompDistinctMZWeak:
        case PartitionType::CompDistinctZero: {
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

int PartitionsCount(const std::vector<int> &Reps, PartDesign &part, int lenV) {

    part.count = 0.0;
    part.numUnknown = false;
    part.bigCount = 0;

    if (part.ptype == PartitionType::NoSolution) {
        return 1;
    }

    const int strtLen = std::count_if(
        part.startZ.cbegin(), part.startZ.cend(), [](int i){return i > 0;}
    );

    const auto no_algo_it = std::find(
        NoCountAlgoPTypeArr.cbegin(), NoCountAlgoPTypeArr.cend(), part.ptype
    );

    if (no_algo_it != NoCountAlgoPTypeArr.end()) {
        part.numUnknown = true;
        return 0;
    }

    if (part.ptype == PartitionType::LengthOne) {
        part.count = static_cast<int>(part.solnExist);
        return 1;
    }

    std::unique_ptr<CountClass> Counter = MakeCount(part.ptype);

    if (!Counter) {
        Counter = MakeMultisetCount(part.ptype, Reps);
    }

    if (Counter) {
        std::vector<int> allowed(part.cap);
        std::iota(allowed.begin(), allowed.end(), 1);

        part.count = Counter->GetCount(
            part.mapTar, part.width, allowed, strtLen
        );

        if (IsBeyondBound(part.count)) {
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
