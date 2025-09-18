#include "Partitions/PartitionsCountMultiset.h"
#include "Partitions/PartitionsCountDistinct.h"
#include "Partitions/BigPartsCountDistinct.h"
#include "Partitions/PartitionsCountRep.h"
#include "Partitions/BigPartsCountRep.h"
#include "Partitions/PartitionsTypes.h"
#include "Partitions/PartitionsCount.h"
#include <numeric>  // std::accumulate

#include "cpp11/R.hpp"
#include "cpp11/protect.hpp"

// The variable k is strtLen
using nthPartsPtr = std::vector<int> (*const)(int n, int m, int cap, int k,
                                      double dblIdx, const mpz_class &mpzIdx);

//*********************** Compositions Functions **************************//

std::vector<int> nthCompsRep(int n, int m, int cap, int k,
                             double dblIdx, const mpz_class &mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    --n;
    --m;

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m, j = 0) {
        for (double temp = CountCompsRepLen(n, m); temp <= dblIdx; ++j) {
            --n;
            dblIdx -= temp;
            temp = CountCompsRepLen(n, m);
        }

        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

std::vector<int> nthCompsRepZero(int n, int m, int cap, int k,
                                 double dblIdx, const mpz_class &mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    bool incr_j = false;
    --m;

    for (int i = 0, j = 0; i < (width - 1); ++i, --m, j = incr_j) {
        double temp = incr_j ? CountCompsRepLen(n, m) :
                    CountCompsRepZNotWk(n, m);

        for (; temp <= dblIdx; ++j) {
            incr_j = true;
            --n;
            dblIdx -= temp;
            temp = CountCompsRepLen(n, m);
        }

        if (incr_j) --n;
        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), 0);
    return res;
}

std::vector<int> nthCompsDistinct(int n, int m, int cap, int k,
                                  double dblIdx, const mpz_class &mpzIdx) {

    const int width = m;
    const int max_val = n - (width * (width - 1)) / 2;

    std::vector<char> mask(n + 1, 0);
    std::vector<int> res(width, 0);
    --m;

    mask[1] = 1;
    std::vector<int> allowed(max_val - 1);
    std::iota(allowed.begin(), allowed.end(), 2);

    int partial_sum = 1;
    int cur_val = 1;

    for (int i = 0, j = 0; i < (width - 1); ++i, --m) {
        double temp = CountCompDistLenRstrctd(n - partial_sum, m, allowed);

        for (; temp <= dblIdx; cur_val = j) {
            while (mask[j + 1]) {
                ++j;
            }

            partial_sum += (j + 1 - cur_val);
            UpdateAllowed(mask, allowed, i, j + 1, width,
                          n, cur_val, partial_sum);

            dblIdx -= temp;
            temp = CountCompDistLenRstrctd(n - partial_sum, m, allowed);

            if (temp <= dblIdx) {
                ++j;
            }
        }

        res[i] = j;
        j = 0;

        while (mask[j + 1]) {
            ++j;
        }

        cur_val = j + 1;
        partial_sum += (j + 1);
        UpdateAllowed(mask, allowed, i + 1, j + 1, width,
                      n, cur_val, partial_sum);
    }

    res[width - 1] = n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

//************************* Partition Functions ***************************//

std::vector<int> nthPartsRepLen(int n, int m, int cap, int k,
                                double dblIdx, const mpz_class &mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    --n;
    --m;

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m) {
        for (double temp = CountPartsRepLen(n, m); temp <= dblIdx; ++j) {
            n -= (m + 1);
            dblIdx -= temp;
            temp = CountPartsRepLen(n, m);
        }

        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

std::vector<int> nthPartsRepShort(int n, int m, int cap, int k,
                                  double dblIdx, const mpz_class &mpzIdx) {

    return nthPartsRepLen(n, m, cap, k, dblIdx, mpzIdx);
}

std::vector<int> nthPartsRep(int n, int m, int cap, int k,
                             double dblIdx, const mpz_class &mpzIdx) {

    return nthPartsRepLen(n * 2, m, cap, k, dblIdx, mpzIdx);
}

std::vector<int> nthPartsRepCap(int n, int m, int cap, int k,
                                double dblIdx, const mpz_class &mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    --n;
    --m;

    std::vector<int> allowed(cap);
    std::iota(allowed.begin(), allowed.end(), 1);

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m) {
        for (double temp = CountPartsRepLenRstrctd(n, m, allowed);
             temp <= dblIdx; ++j) {
            n -= (m + 1);
            allowed.pop_back();
            dblIdx -= temp;
            temp = CountPartsRepLenRstrctd(n, m, allowed);
        }

        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

std::vector<int> nthPartsDistinctLen(int n, int m, int cap, int k,
                                     double dblIdx, const mpz_class &mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    n -= m;
    --m;

    for (int i = 0, j = 0; i < (width - 1); ++i, n -= m, --m, ++j) {
        for (double temp = CountPartsDistinctLen(n, m); temp <= dblIdx; ++j) {
            n -= (m + 1);
            dblIdx -= temp;
            temp = CountPartsDistinctLen(n, m);
        }

        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

std::vector<int> nthPartsDistinctOneZero(
    int n, int m, int cap, int k, double dblIdx, const mpz_class &mpzIdx
) {

    return nthPartsDistinctLen(n, m, cap, k, dblIdx, mpzIdx);
}

std::vector<int> nthPartsDistinctMultiZero(
    int n, int m, int cap, int k, double dblIdx, const mpz_class &mpzIdx
) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    bool incr_j = false;
    --m;

    std::vector<int> empty_allowed;

    for (int i = 0, j = 0; i < (width - 1); ++i, --m) {
        double temp = (incr_j || i >= (width - k)) ?
                      CountPartsDistinctLen(n, m) :
                      CountPartsDistinctMultiZero(n, m, empty_allowed, k);

        for (; temp <= dblIdx; ++j) {
            incr_j = true;
            n -= (m + 1);
            dblIdx -= temp;
            temp = CountPartsDistinctLen(n, m);
        }

        res[i] = j;

        if (incr_j || (i + 1) >= (width - k)) {
            ++j;
            n -= m;
        }
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), 0);
    return res;
}

std::vector<int> nthPartsDistinctCap(int n, int m, int cap, int k,
                                     double dblIdx, const mpz_class &mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    n -= m;
    --cap;
    --m;

    std::vector<int> allowed(cap);
    std::iota(allowed.begin(), allowed.end(), 1);

    for (int i = 0, j = 0; i < (width - 1); ++i, n -= m, --m, ++j) {
        for (double temp = CountPartsDistLenRstrctd(n, m, allowed);
             temp <= dblIdx; ++j) {
            n -= (m + 1);
            allowed.pop_back();
            dblIdx -= temp;
            temp = CountPartsDistLenRstrctd(n, m, allowed);
        }

        res[i] = j;
        allowed.pop_back();
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

std::vector<int> nthPartsDistinctCapMZ(
    int n, int m, int cap, int k, double dblIdx, const mpz_class &mpzIdx
) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    bool incr_j = false;
    --m;

    std::vector<int> allowed(cap);
    std::iota(allowed.begin(), allowed.end(), 1);

    for (int i = 0, j = 0; i < (width - 1); ++i, --m) {
        double temp = (incr_j || i >= (width - k)) ?
                CountPartsDistLenRstrctd(n, m, allowed) :
                CountPartsDistinctRstrctdMZ(n, m, allowed, k);

        for (; temp <= dblIdx; ++j) {
            incr_j = true;
            n -= (m + 1);
            allowed.pop_back();
            dblIdx -= temp;
            temp = CountPartsDistLenRstrctd(n, m, allowed);
        }

        res[i] = j;

        if (incr_j || (i + 1) >= (width - k)) {
            ++j;
            n -= m;
            allowed.pop_back();
        }
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), 0);
    return res;
}

//*********************** Starting Gmp Functions **************************//

std::vector<int> nthCompsRepGmp(int n, int m, int cap, int k,
                                double dblIdx, const mpz_class &mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    --n;
    --m;

    mpz_class temp;
    mpz_class index(mpzIdx);

    std::unique_ptr<CountClass> Counter = MakeCount(
        PartitionType::CompRepNoZero
    );

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m, j = 0) {
        Counter->GetCount(temp, n, m);

        for (; cmp(temp, index) <= 0; ++j) {
            --n;
            index -= temp;
            Counter->GetCount(temp, n, m);
        }

        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

std::vector<int> nthCompsRepZeroGmp(int n, int m, int cap, int k,
                                    double dblIdx, const mpz_class &mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    bool incr_j = false;
    --m;

    mpz_class temp;
    mpz_class index(mpzIdx);
    std::unique_ptr<CountClass> Counter = MakeCount(
        PartitionType::CmpRpZroNotWk
    );

    std::vector<int> empty_allowed;

    for (int i = 0, j = 0; i < (width - 1); ++i, --m, j = incr_j) {
        Counter->GetCount(temp, n, m, empty_allowed, k, !incr_j);

        for (; cmp(temp, index) <= 0; ++j) {
            incr_j = true;
            --n;
            index -= temp;
            Counter->GetCount(temp, n, m, empty_allowed, k, false);
        }

        temp = 0;
        if (incr_j) --n;
        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), 0);
    return res;
}

std::vector<int> nthCompsDistinctGmp(int n, int m, int cap, int k,
                                     double dblIdx, const mpz_class &mpzIdx) {

    const int width = m;
    const int max_val = n - (width * (width - 1)) / 2;

    std::vector<char> mask(n + 1, 0);
    std::vector<int> res(width, 0);
    --m;

    mask[1] = 1;
    std::vector<int> allowed(max_val - 1);
    std::iota(allowed.begin(), allowed.end(), 2);

    mpz_class temp;
    mpz_class index(mpzIdx);

    const PartitionType ptype = PartitionType::PrmDstPrtCap;
    std::unique_ptr<CountClass> Counter = MakeCount(ptype);

    Counter->SetArrSize(ptype, n, m);
    Counter->InitializeMpz();

    int partial_sum = 1;
    int cur_val = 1;

    for (int i = 0, j = 0; i < (width - 1); ++i, --m) {
        Counter->GetCount(temp, n - partial_sum, m, allowed);

        for (; cmp(temp, index) <= 0; cur_val = j) {
            while (mask[j + 1]) {
                ++j;
            }

            partial_sum += (j + 1 - cur_val);
            UpdateAllowed(mask, allowed, i, j + 1, width,
                          n, cur_val, partial_sum);

            index -= temp;
            Counter->GetCount(temp, n - partial_sum, m, allowed);

            if (cmp(temp, index) <= 0) {
                ++j;
            }
        }

        res[i] = j;
        j = 0;

        while (mask[j + 1]) {
            ++j;
        }

        cur_val = j + 1;
        partial_sum += (j + 1);

        UpdateAllowed(mask, allowed, i + 1, j + 1, width,
                      n, cur_val, partial_sum);
    }

    res[width - 1] = n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

std::vector<int> nthPartsRepLenGmp(int n, int m, int cap, int k,
                                   double dblIdx, const mpz_class &mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    --n;
    --m;

    mpz_class temp;
    mpz_class index(mpzIdx);

    const PartitionType ptype = PartitionType::RepShort;
    std::unique_ptr<CountClass> Counter = MakeCount(ptype);

    Counter->SetArrSize(ptype, n, m);
    Counter->InitializeMpz();

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m) {
        Counter->GetCount(temp, n, m);

        for (; cmp(temp, index) <= 0; ++j) {
            n -= (m + 1);
            index -= temp;
            Counter->GetCount(temp, n, m);
        }

        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

std::vector<int> nthPartsRepShortGmp(int n, int m, int cap, int k,
                                     double dblIdx, const mpz_class &mpzIdx) {

    return nthPartsRepLenGmp(n, m, cap, k, dblIdx, mpzIdx);
}

std::vector<int> nthPartsRepGmp(int n, int m, int cap, int k,
                                double dblIdx, const mpz_class &mpzIdx) {

    return nthPartsRepLenGmp(n * 2, m, cap, k, dblIdx, mpzIdx);
}

std::vector<int> nthPartsRepCapGmp(int n, int m, int cap, int k,
                                   double dblIdx, const mpz_class &mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    --n;
    --m;

    mpz_class temp;
    mpz_class index(mpzIdx);

    const PartitionType ptype = PartitionType::RepCapped;
    std::unique_ptr<CountClass> Counter = MakeCount(ptype);

    Counter->SetArrSize(ptype, n, m);
    Counter->InitializeMpz();

    std::vector<int> allowed(cap);
    std::iota(allowed.begin(), allowed.end(), 1);

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m) {
        Counter->GetCount(temp, n, m, allowed, k);

        for (; cmp(temp, index) <= 0; ++j) {
            n -= (m + 1);
            allowed.pop_back();
            index -= temp;
            Counter->GetCount(temp, n, m, allowed, k);
        }

        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

std::vector<int> nthPartsDistinctLenGmp(
    int n, int m, int cap, int k, double dblIdx, const mpz_class &mpzIdx
) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    n -= m;
    --m;

    mpz_class temp;
    mpz_class index(mpzIdx);

    const PartitionType ptype = PartitionType::DstctNoZero;
    std::unique_ptr<CountClass> Counter = MakeCount(ptype);

    Counter->SetArrSize(ptype, n, m);
    Counter->InitializeMpz();

    for (int i = 0, j = 0; i < (width - 1); ++i, n -= m, --m, ++j) {
        Counter->GetCount(temp, n, m);

        for (; cmp(temp, index) <= 0; ++j) {
            n -= (m + 1);
            index -= temp;
            Counter->GetCount(temp, n, m);
        }

        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

std::vector<int> nthPartsDistinctOneZeroGmp(
    int n, int m, int cap, int k, double dblIdx, const mpz_class &mpzIdx
) {

    return nthPartsDistinctLenGmp(n, m, cap, k, dblIdx, mpzIdx);
}

std::vector<int> nthPartsDistinctMultiZeroGmp(
    int n, int m, int cap, int k, double dblIdx, const mpz_class &mpzIdx
) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    bool incr_j = false;
    --m;

    mpz_class temp;
    mpz_class index(mpzIdx);

    const PartitionType ptype = PartitionType::DstctMultiZero;
    std::unique_ptr<CountClass> Counter = MakeCount(ptype);

    Counter->SetArrSize(ptype, n, m);
    Counter->InitializeMpz();

    std::vector<int> empty_allowed;

    for (int i = 0, j = 0; i < (width - 1); ++i, --m) {
        const bool bLiteral = !(incr_j || i >= (width - k));
        Counter->GetCount(temp, n, m, empty_allowed, k, bLiteral);

        for (; cmp(temp, index) <= 0; ++j) {
            incr_j = true;
            n -= (m + 1);
            index -= temp;
            Counter->GetCount(temp, n, m, empty_allowed, k, false);
        }

        res[i] = j;

        if (incr_j || (i + 1) >= (width - k)) {
            ++j;
            n -= m;
        }
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), 0);
    return res;
}

std::vector<int> nthPartsDistinctCapGmp(
    int n, int m, int cap, int k, double dblIdx, const mpz_class &mpzIdx
) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    n -= m;
    --cap;
    --m;

    mpz_class temp;
    mpz_class index(mpzIdx);

    const PartitionType ptype = PartitionType::DstctCapped;
    std::unique_ptr<CountClass> Counter = MakeCount(ptype);

    Counter->SetArrSize(ptype, n, m);
    Counter->InitializeMpz();

    std::vector<int> allowed(cap);
    std::iota(allowed.begin(), allowed.end(), 1);

    for (int i = 0, j = 0; i < (width - 1); ++i, n -= m, --m, ++j) {
        Counter->GetCount(temp, n, m, allowed, k);

        for (; cmp(temp, index) <= 0; ++j) {
            n -= (m + 1);
            allowed.pop_back();
            index -= temp;
            Counter->GetCount(temp, n, m, allowed, k);
        }

        res[i] = j;
        allowed.pop_back();
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

std::vector<int> nthPartsDistinctCapMZGmp(
    int n, int m, int cap, int k, double dblIdx, const mpz_class &mpzIdx
) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    bool incr_j = false;
    --m;

    mpz_class temp;
    mpz_class index(mpzIdx);

    const PartitionType ptype = PartitionType::DstctCappedMZ;
    std::unique_ptr<CountClass> Counter = MakeCount(ptype);

    Counter->SetArrSize(ptype, n, m);
    Counter->InitializeMpz();

    std::vector<int> allowed(cap);
    std::iota(allowed.begin(), allowed.end(), 1);

    for (int i = 0, j = 0; i < (width - 1); ++i, --m) {
        const bool bLiteral = !(incr_j || i >= (width - k));
        Counter->GetCount(temp, n, m, allowed, k, bLiteral);

        for (; cmp(temp, index) <= 0; ++j) {
            incr_j = true;
            n -= (m + 1);
            allowed.pop_back();
            index -= temp;
            Counter->GetCount(temp, n, m, allowed, k, false);
        }

        res[i] = j;

        if (incr_j || (i + 1) >= (width - k)) {
            ++j;
            n -= m;
            allowed.pop_back();
        }
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), 0);
    return res;
}

nthPartsPtr GetNthPartsFunc(PartitionType ptype, bool IsGmp) {

    if (IsGmp) {
        switch (ptype) {
            case PartitionType::DstctCapped: {
                return(nthPartsPtr(nthPartsDistinctCapGmp));
            } case PartitionType::DstctCappedMZ: {
                return(nthPartsPtr(nthPartsDistinctCapMZGmp));
            } case PartitionType::DstctNoZero: {
                return(nthPartsPtr(nthPartsDistinctLenGmp));
            } case PartitionType::DstctOneZero: {
                return(nthPartsPtr(nthPartsDistinctOneZeroGmp));
            } case PartitionType::DstctMultiZero: {
                return(nthPartsPtr(nthPartsDistinctMultiZeroGmp));
            } case PartitionType::DstctStdAll: {
                return(nthPartsPtr(nthPartsDistinctMultiZeroGmp));
            } case PartitionType::RepCapped: {
                return(nthPartsPtr(nthPartsRepCapGmp));
            } case PartitionType::RepNoZero: {
                return(nthPartsPtr(nthPartsRepLenGmp));
            } case PartitionType::RepShort: {
                return(nthPartsPtr(nthPartsRepShortGmp));
            } case PartitionType::RepStdAll: {
                return(nthPartsPtr(nthPartsRepGmp));
            } case PartitionType::CompRepNoZero: {
                return(nthPartsPtr(nthCompsRepGmp));
            } case PartitionType::CompRepWeak: {
                return(nthPartsPtr(nthCompsRepGmp));
            } case PartitionType::CmpRpZroNotWk: {
                return(nthPartsPtr(nthCompsRepZeroGmp));
            } case PartitionType::CmpDstctNoZero: {
                return(nthPartsPtr(nthCompsDistinctGmp));
            } default : {
                cpp11::stop("No algorithm available");
            }
        }
    } else {
        switch (ptype) {
            case PartitionType::DstctCapped: {
                return(nthPartsPtr(nthPartsDistinctCap));
            } case PartitionType::DstctCappedMZ: {
                return(nthPartsPtr(nthPartsDistinctCapMZ));
            } case PartitionType::DstctNoZero: {
                return(nthPartsPtr(nthPartsDistinctLen));
            } case PartitionType::DstctOneZero: {
                return(nthPartsPtr(nthPartsDistinctOneZero));
            } case PartitionType::DstctMultiZero: {
                return(nthPartsPtr(nthPartsDistinctMultiZero));
            } case PartitionType::DstctStdAll: {
                return(nthPartsPtr(nthPartsDistinctMultiZero));
            } case PartitionType::RepCapped: {
                return(nthPartsPtr(nthPartsRepCap));
            } case PartitionType::RepNoZero: {
                return(nthPartsPtr(nthPartsRepLen));
            } case PartitionType::RepShort: {
                return(nthPartsPtr(nthPartsRepShort));
            } case PartitionType::RepStdAll: {
                return(nthPartsPtr(nthPartsRep));
            } case PartitionType::CompRepNoZero: {
                return(nthPartsPtr(nthCompsRep));
            } case PartitionType::CompRepWeak: {
                return(nthPartsPtr(nthCompsRep));
            } case PartitionType::CmpRpZroNotWk: {
                return(nthPartsPtr(nthCompsRepZero));
            } case PartitionType::CmpDstctNoZero: {
                return(nthPartsPtr(nthCompsDistinct));
            } default : {
                cpp11::stop("No algorithm available");
            }
        }
    }
}
