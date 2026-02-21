#include "cpp11/R.hpp"
#include "cpp11/protect.hpp"

#include "Partitions/PartitionsCountMultiset.h"
#include "Partitions/PartitionsCountDistinct.h"
#include "Partitions/BigPartsCountDistinct.h"
#include "Partitions/PartitionsCountRep.h"
#include "Partitions/BigPartsCountRep.h"
#include "Partitions/PartitionsCount.h"
#include "Partitions/RankPartition.h"
#include "CppConvert/Constants.h"
#include <numeric>  // std::accumulate

//*********************** Trivial Length One Case **************************//

void rankLengthOne(std::vector<int>::iterator iter, int n, int m,
                   int cap, int k, double &dblIdx, mpz_class &mpzIdx) {
    dblIdx = 0;
}

//*********************** Compositions Functions **************************//

void rankCompsRep(std::vector<int>::iterator iter, int n, int m,
                  int cap, int k, double &dblIdx, mpz_class &mpzIdx) {

    const int width = m;
    dblIdx = 0;

    --n;
    --m;

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m, j = 0, ++iter) {
        double temp = CountCompsRepLen(n, m);

        for (int idx = *iter; j < idx; ++j) {
            --n;
            dblIdx += temp;
            temp = CountCompsRepLen(n, m);
        }
    }
}

void rankCompsRepZero(std::vector<int>::iterator iter, int n, int m,
                      int cap, int k, double &dblIdx, mpz_class &mpzIdx) {

    const int width = m;
    dblIdx = 0;

    bool incr_j = false;
    --m;

    for (int i = 0, j = 0; i < (width - 1); ++i, --m, j = incr_j, ++iter) {
        double temp = incr_j ? CountCompsRepLen(n, m) :
            CountCompsRepZNotWk(n, m);

        for (int idx = *iter; j < idx; ++j) {
            incr_j = true;
            --n;
            dblIdx += temp;
            temp = CountCompsRepLen(n, m);
        }

        if (incr_j) --n;
    }
}

void rankCompsDistinct(std::vector<int>::iterator iter, int n, int m,
                       int cap, int k, double &dblIdx, mpz_class &mpzIdx) {

    const int width = m;
    const int max_val = std::min(cap, n - (width * (width - 1)) / 2);

    std::vector<char> mask(max_val + 1, 0);
    const int mask_size = mask.size();
    --m;

    mask[1] = 1;
    std::vector<int> allowed(max_val - 1);
    std::iota(allowed.begin(), allowed.end(), 2);

    int partial_sum = 1;
    int cur_val = 1;

    for (int i = 0, j = 0; i < (width - 1); ++i, --m, ++iter) {
        double temp = CountCompDistLenRstrctd(n - partial_sum, m, allowed);
        bool keepGoing = false;

        for (int idx = *iter; j < idx || keepGoing; cur_val = j) {
            while ((j + 1) < mask_size && mask[j + 1]) {
                ++j;
            }

            partial_sum += (j + 1 - cur_val);
            UpdateAllowed(mask, allowed, i, j + 1, width,
                          n, cur_val, partial_sum);

            dblIdx += temp;
            temp = CountCompDistLenRstrctd(n - partial_sum, m, allowed);

            ++j;
            keepGoing = j <= idx;
        }

        j = 0;

        while ((j + 1) < mask_size && mask[j + 1]) {
            ++j;
        }

        cur_val = j + 1;
        partial_sum += (j + 1);
        UpdateAllowed(mask, allowed, i + 1, j + 1, width,
                      n, cur_val, partial_sum);
    }
}

void rankCompsDistinctMZ(std::vector<int>::iterator iter, int n, int m,
                         int cap, int k, double &dblIdx, mpz_class &mpzIdx) {

    std::vector<int> allowed(cap);
    std::iota(allowed.begin(), allowed.end(), 1);
    dblIdx = 0;

    double temp = (cap == n) ?
        CountCompsDistinctLen(n, k) :
            CountCompDistLenRstrctd(n, k, allowed);

    // NOTE: iter is guaranteed by the calling rank dispatcher to reference
    // a contiguous buffer of at least m elements representing the index vector.
    // Therefore [iter, iter + m) is always a valid range.
    std::vector<int> z(iter, iter + m);  // snapshot
    for (int &x : z) --x;                 // shift encoding

    auto z_it = z.begin();
    int width = m;

    for (int zz = k; *z_it == -1 && zz < m; ++zz, ++z_it) {
        --width;
    }

    while (k < width) {
        dblIdx += temp;
        ++k;
        temp = (cap == n) ?
            CountCompsDistinctLen(n, k) :
                CountCompDistLenRstrctd(n, k, allowed);
    }

    rankCompsDistinct(z_it, n, k, cap, k, dblIdx, mpzIdx);
}

void rankCompsDistinctMZWeak(
    std::vector<int>::iterator iter, int n, int m, int cap, int k,
    double &dblIdx, mpz_class &mpzIdx
) {

    const int width = m;
    const int max_val = std::min(cap, n - (k * (k - 1)) / 2);
    int zeros_remaining = width - k;

    std::vector<char> mask(max_val + 1, 0);
    --zeros_remaining;
    --m;

    mask[0] = 1;
    std::vector<int> allowed(max_val);
    std::iota(allowed.begin(), allowed.end(), 1);

    int partial_sum = 0;
    int cur_val = 0;

    bool anyZeros = true;
    const int mask_size = mask.size();

    for (int i = 0, j = 0; i < (width - 1); ++i, --m, ++iter) {
        int strtLen = std::max(1, m - zeros_remaining);
        if (anyZeros) ++zeros_remaining;

        double temp = CountCompsDistinctRstrctdMZWeak(
            n - partial_sum, m, allowed, strtLen
        );

        for (int idx = *iter; j < idx && n > partial_sum; cur_val = j) {
            ++j;

            while (j < mask_size && mask[j]) {
                ++j;
            }

            partial_sum += (j - cur_val);
            UpdateAllowed(mask, allowed, i, j, width - zeros_remaining,
                          n, cur_val, partial_sum);

            dblIdx += temp;
            strtLen = std::max(1, m - zeros_remaining);
            temp = CountCompsDistinctRstrctdMZWeak(
                n - partial_sum, m, allowed, strtLen
            );
        }

        if (j == 0) --zeros_remaining;
        anyZeros = zeros_remaining > 0;
        j = anyZeros ? 0 : 1;

        while (!anyZeros && j < mask_size && mask[j]) {
            ++j;
        }

        if (anyZeros) --zeros_remaining;
        cur_val = j;
        partial_sum += j;
        UpdateAllowed(mask, allowed, i + 1, j, width - zeros_remaining,
                      n, cur_val, partial_sum);
    }
}

//************************* Partition Functions ***************************//

void rankPartsRepLen(std::vector<int>::iterator iter, int n, int m,
                     int cap, int k, double &dblIdx, mpz_class &mpzIdx) {

    const int width = m;
    dblIdx = 0;

    --n;
    --m;

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m, ++iter) {
        double temp = CountPartsRepLen(n, m);

        for (int idx = *iter; j < idx; ++j) {
            n -= (m + 1);
            dblIdx += temp;
            temp = CountPartsRepLen(n, m);
        }
    }
}

void rankPartsRepShort(std::vector<int>::iterator iter, int n, int m,
                       int cap, int k, double &dblIdx, mpz_class &mpzIdx) {

    rankPartsRepLen(iter, n, m, cap, k, dblIdx, mpzIdx);
}

void rankPartsRep(std::vector<int>::iterator iter, int n, int m,
                  int cap, int k, double &dblIdx, mpz_class &mpzIdx) {

    rankPartsRepLen(iter, n * 2, m, cap, k, dblIdx, mpzIdx);
}

void rankPartsRepCap(std::vector<int>::iterator iter, int n, int m,
                     int cap, int k, double &dblIdx, mpz_class &mpzIdx) {

    const int width = m;
    dblIdx = 0;

    --n;
    --m;

    std::vector<int> allowed(cap);
    std::iota(allowed.begin(), allowed.end(), 1);

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m, ++iter) {
        double temp = CountPartsRepLenRstrctd(n, m, allowed);

        for (int idx = *iter; j < idx; ++j) {
            n -= (m + 1);
            allowed.pop_back();
            dblIdx += temp;
            temp = CountPartsRepLenRstrctd(n, m, allowed);
        }
    }
}

void rankPartsDistinctLen(std::vector<int>::iterator iter, int n, int m,
                          int cap, int k, double &dblIdx, mpz_class &mpzIdx) {

    const int width = m;
    dblIdx = 0;

    n -= m;
    --m;

    for (int i = 0, j = 0; i < (width - 1); ++i, n -= m, --m, ++j, ++iter) {
        double temp = CountPartsDistinctLen(n, m);

        for (int idx = *iter; j < idx; ++j) {
            n -= (m + 1);
            dblIdx += temp;
            temp = CountPartsDistinctLen(n, m);
        }
    }
}

void rankPartsDistinctOneZero(
    std::vector<int>::iterator iter, int n, int m,
    int cap, int k, double &dblIdx, mpz_class &mpzIdx
) {

    rankPartsDistinctLen(iter, n, m, cap, k, dblIdx, mpzIdx);
}

void rankPartsDistinctMultiZero(
    std::vector<int>::iterator iter, int n, int m,
    int cap, int k, double &dblIdx, mpz_class &mpzIdx
) {

    const int width = m;
    dblIdx = 0;

    bool incr_j = false;
    --m;

    std::vector<int> empty_allowed;

    for (int i = 0, j = 0; i < (width - 1); ++i, --m, ++iter) {
        double temp = (incr_j || i >= (width - k)) ?
                      CountPartsDistinctLen(n, m) :
                      CountPartsDistinctMultiZero(n, m, empty_allowed, k);

        for (int idx = *iter; j < idx; ++j) {
            incr_j = true;
            n -= (m + 1);
            dblIdx += temp;
            temp = CountPartsDistinctLen(n, m);
        }

        if (incr_j || (i + 1) >= (width - k)) {
            ++j;
            n -= m;
        }
    }
}

void rankPartsDistinctCap(std::vector<int>::iterator iter, int n, int m,
                          int cap, int k, double &dblIdx, mpz_class &mpzIdx) {

    const int width = m;
    dblIdx = 0;

    n -= m;
    --cap;
    --m;

    std::vector<int> allowed(cap);
    std::iota(allowed.begin(), allowed.end(), 1);

    for (int i = 0, j = 0; i < (width - 1); ++i, n -= m, --m, ++j, ++iter) {
        double temp = CountPartsDistLenRstrctd(n, m, allowed, k);

        for (int idx = *iter; j < idx; ++j) {
            n -= (m + 1);
            allowed.pop_back();
            dblIdx += temp;
            temp = CountPartsDistLenRstrctd(n, m, allowed, k);
        }

        allowed.pop_back();
    }
}

void rankPartsDistinctCapMZ(std::vector<int>::iterator iter, int n, int m,
                            int cap, int k, double &dblIdx, mpz_class &mpzIdx) {

    const int width = m;
    dblIdx = 0;

    bool incr_j = false;
    --m;

    std::vector<int> allowed(cap);
    std::iota(allowed.begin(), allowed.end(), 1);

    for (int i = 0, j = 0; i < (width - 1); ++i, --m, ++iter) {
        double temp = (incr_j || i >= (width - k)) ?
                    CountPartsDistLenRstrctd(n, m, allowed, k) :
                    CountPartsDistinctRstrctdMZ(n, m, allowed, k);

        for (int idx = *iter; j < idx; ++j) {
            incr_j = true;
            n -= (m + 1);
            allowed.pop_back();
            dblIdx += temp;
            temp = CountPartsDistLenRstrctd(n, m, allowed, k);
        }

        if (incr_j || (i + 1) >= (width - k)) {
            ++j;
            n -= m;
            allowed.pop_back();
        }
    }
}

//*********************** Starting Gmp Functions **************************//

void rankCompsRepGmp(std::vector<int>::iterator iter, int n, int m,
                     int cap, int k, double &dblIdx, mpz_class &mpzIdx) {

    const int width = m;
    mpz_class temp;
    mpzIdx = 0;

    --n;
    --m;

    std::unique_ptr<CountClass> Counter = MakeCount(
        PartitionType::CompRepNoZero
    );

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m, j = 0, ++iter) {
        Counter->GetCount(temp, n, m);

        for (int idx = *iter; j < idx; ++j) {
            --n;
            mpzIdx += temp;
            Counter->GetCount(temp, n, m);
        }
    }
}

void rankCompsRepZeroGmp(std::vector<int>::iterator iter, int n, int m,
                         int cap, int k, double &dblIdx, mpz_class &mpzIdx) {

    const int width = m;
    mpzIdx = 0;

    bool incr_j = false;
    --m;

    mpz_class temp;
    std::unique_ptr<CountClass> Counter = MakeCount(
        PartitionType::CmpRpZroNotWk
    );

    std::vector<int> empty_allowed;

    for (int i = 0, j = 0; i < (width - 1); ++i, --m, j = incr_j, ++iter) {
        Counter->GetCount(temp, n, m, empty_allowed, k, !incr_j);

        for (int idx = *iter; j < idx; ++j) {
            incr_j = true;
            --n;
            mpzIdx += temp;
            Counter->GetCount(temp, n, m, empty_allowed, k, false);
        }

        temp = 0;
        if (incr_j) --n;
    }
}

void rankCompsDistinctGmp(std::vector<int>::iterator iter, int n, int m,
                          int cap, int k, double &dblIdx, mpz_class &mpzIdx) {

    const int width = m;
    const int max_val = std::min(cap, n - (width * (width - 1)) / 2);

    std::vector<char> mask(max_val + 1, 0);
    const int mask_size = mask.size();
    --m;

    mask[1] = 1;
    std::vector<int> allowed(max_val - 1);
    std::iota(allowed.begin(), allowed.end(), 2);

    mpz_class temp;
    const PartitionType ptype = PartitionType::PrmDstPrtCap;
    std::unique_ptr<CountClass> Counter = MakeCount(ptype);

    Counter->SetArrSize(ptype, n, m);
    Counter->InitializeMpz();

    int partial_sum = 1;
    int cur_val = 1;

    for (int i = 0, j = 0; i < (width - 1); ++i, --m, ++iter) {
        Counter->GetCount(temp, n - partial_sum, m, allowed);
        bool keepGoing = false;

        for (int idx = *iter; j < idx || keepGoing; cur_val = j) {
            while ((j + 1) < mask_size && mask[j + 1]) {
                ++j;
            }

            partial_sum += (j + 1 - cur_val);
            UpdateAllowed(mask, allowed, i, j + 1, width,
                          n, cur_val, partial_sum);

            mpzIdx += temp;
            Counter->GetCount(temp, n - partial_sum, m, allowed);

            ++j;
            keepGoing = j <= idx;
        }

        j = 0;

        while ((j + 1) < mask_size && mask[j + 1]) {
            ++j;
        }

        cur_val = j + 1;
        partial_sum += (j + 1);
        UpdateAllowed(mask, allowed, i + 1, j + 1, width,
                      n, cur_val, partial_sum);
    }
}

void rankCompsDistinctMZGmp(std::vector<int>::iterator iter, int n, int m,
                            int cap, int k, double &dblIdx, mpz_class &mpzIdx) {

    mpz_class temp;

    // NOTE: iter is guaranteed by the calling rank dispatcher to reference
    // a contiguous buffer of at least m elements representing the index vector.
    // Therefore [iter, iter + m) is always a valid range.
    std::vector<int> z(iter, iter + m);  // snapshot
    for (int &x : z) --x;                 // shift encoding

    auto z_it = z.begin();
    int width = m;

    for (int zz = k; *z_it == -1 && zz < m; ++zz, ++z_it) {
        --width;
    }

    if (cap == n) {
        const PartitionType ptype = PartitionType::CmpDstctNoZero;
        std::unique_ptr<CountClass> Counter = MakeCount(ptype);

        Counter->SetArrSize(ptype, n, m);
        Counter->InitializeMpz();
        Counter->GetCount(temp, n, k);

        while (k < width) {
            mpzIdx += temp;
            ++k;
            Counter->GetCount(temp, n, k);
        }
    } else {
        std::vector<int> allowed(cap);
        std::iota(allowed.begin(), allowed.end(), 1);

        const PartitionType ptype = PartitionType::CmpDstctCapped;
        std::unique_ptr<CountClass> Counter = MakeCount(ptype);

        Counter->SetArrSize(ptype, n, m);
        Counter->InitializeMpz();
        Counter->GetCount(temp, n, k, allowed);

        while (k < width) {
            mpzIdx += temp;
            ++k;
            Counter->GetCount(temp, n, k, allowed);
        }
    }

    rankCompsDistinctGmp(z_it, n, k, cap, k, dblIdx, mpzIdx);
}

void rankCompsDistinctMZWeakGmp(
    std::vector<int>::iterator iter, int n, int m, int cap, int k,
    double &dblIdx, mpz_class &mpzIdx
) {

    const int width = m;
    const int max_val = std::min(cap, n - (k * (k - 1)) / 2);
    int zeros_remaining = width - k;

    std::vector<char> mask(max_val + 1, 0);
    --zeros_remaining;
    --m;

    mask[0] = 1;
    std::vector<int> allowed(max_val);
    std::iota(allowed.begin(), allowed.end(), 1);

    int partial_sum = 0;
    int cur_val = 0;

    bool anyZeros = true;
    const int mask_size = mask.size();

    mpz_class temp;
    const PartitionType ptype = PartitionType::CmpDstCapMZWeak;
    std::unique_ptr<CountClass> Counter = MakeCount(ptype);

    Counter->SetArrSize(ptype, n, m);
    Counter->InitializeMpz();

    for (int i = 0, j = 0; i < (width - 1); ++i, --m, ++iter) {
        int strtLen = std::max(1, m - zeros_remaining);
        if (anyZeros) ++zeros_remaining;

        Counter->GetCount(
            temp, n - partial_sum, m, allowed, strtLen, true
        );

        for (int idx = *iter; j < idx && n > partial_sum; cur_val = j) {
            ++j;

            while (j < mask_size && mask[j]) {
                ++j;
            }

            partial_sum += (j - cur_val);
            UpdateAllowed(mask, allowed, i, j, width - zeros_remaining,
                          n, cur_val, partial_sum);

            mpzIdx += temp;
            strtLen = std::max(1, m - zeros_remaining);
            Counter->GetCount(
                temp, n - partial_sum, m, allowed, strtLen, true
            );
        }

        if (j == 0) --zeros_remaining;
        anyZeros = zeros_remaining > 0;
        j = anyZeros ? 0 : 1;

        while (!anyZeros && j < mask_size && mask[j]) {
            ++j;
        }

        if (anyZeros) --zeros_remaining;
        cur_val = j;
        partial_sum += j;
        UpdateAllowed(mask, allowed, i + 1, j, width - zeros_remaining,
                      n, cur_val, partial_sum);
    }
}

void rankPartsRepLenGmp(std::vector<int>::iterator iter, int n, int m,
                        int cap, int k, double &dblIdx, mpz_class &mpzIdx) {

    const int width = m;
    mpzIdx = 0;

    --n;
    --m;

    mpz_class temp;
    const PartitionType ptype = PartitionType::RepShort;
    std::unique_ptr<CountClass> Counter = MakeCount(ptype);

    Counter->SetArrSize(ptype, n, m);
    Counter->InitializeMpz();

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m, ++iter) {
        Counter->GetCount(temp, n, m);

        for (int idx = *iter; j < idx; ++j) {
            n -= (m + 1);
            mpzIdx += temp;
            Counter->GetCount(temp, n, m);
        }
    }
}

void rankPartsRepShortGmp(std::vector<int>::iterator iter, int n, int m,
                          int cap, int k, double &dblIdx, mpz_class &mpzIdx) {

    rankPartsRepLenGmp(iter, n, m, cap, k, dblIdx, mpzIdx);
}

void rankPartsRepGmp(std::vector<int>::iterator iter, int n, int m,
                     int cap, int k, double &dblIdx, mpz_class &mpzIdx) {

    rankPartsRepLenGmp(iter, n * 2, m, cap, k, dblIdx, mpzIdx);
}

void rankPartsRepCapGmp(std::vector<int>::iterator iter, int n, int m,
                        int cap, int k, double &dblIdx, mpz_class &mpzIdx) {

    const int width = m;
    mpzIdx = 0;

    --n;
    --m;

    mpz_class temp;
    const PartitionType ptype = PartitionType::RepCapped;
    std::unique_ptr<CountClass> Counter = MakeCount(ptype);

    Counter->SetArrSize(ptype, n, m);
    Counter->InitializeMpz();

    std::vector<int> allowed(cap);
    std::iota(allowed.begin(), allowed.end(), 1);

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m, ++iter) {
        Counter->GetCount(temp, n, m, allowed, k);

        for (int idx = *iter; j < idx; ++j) {
            n -= (m + 1);
            allowed.pop_back();
            mpzIdx += temp;
            Counter->GetCount(temp, n, m, allowed, k);
        }
    }
}

void rankPartsDistinctLenGmp(std::vector<int>::iterator iter, int n, int m,
                             int cap, int k, double &dblIdx, mpz_class &mpzIdx) {

    const int width = m;
    mpzIdx = 0;

    n -= m;
    --m;

    mpz_class temp;
    const PartitionType ptype = PartitionType::DstctNoZero;
    std::unique_ptr<CountClass> Counter = MakeCount(ptype);

    Counter->SetArrSize(ptype, n, m);
    Counter->InitializeMpz();

    for (int i = 0, j = 0; i < (width - 1); ++i, n -= m, --m, ++j, ++iter) {
        Counter->GetCount(temp, n, m);

        for (int idx = *iter; j < idx; ++j) {
            n -= (m + 1);
            mpzIdx += temp;
            Counter->GetCount(temp, n, m);
        }
    }
}

void rankPartsDistinctOneZeroGmp(std::vector<int>::iterator iter,
                                 int n, int m, int cap, int k,
                                 double &dblIdx, mpz_class &mpzIdx) {

    rankPartsDistinctLenGmp(iter, n, m, cap, k, dblIdx, mpzIdx);
}

void rankPartsDistinctMultiZeroGmp(std::vector<int>::iterator iter,
                                   int n, int m, int cap, int k,
                                   double &dblIdx, mpz_class &mpzIdx) {

    const int width = m;
    mpzIdx = 0;

    bool incr_j = false;
    --m;

    mpz_class temp;
    const PartitionType ptype = PartitionType::DstctMultiZero;
    std::unique_ptr<CountClass> Counter = MakeCount(ptype);

    Counter->SetArrSize(ptype, n, m);
    Counter->InitializeMpz();

    std::vector<int> empty_allowed;

    for (int i = 0, j = 0; i < (width - 1); ++i, --m, ++iter) {
        const bool bLiteral = !(incr_j || i >= (width - k));
        Counter->GetCount(temp, n, m, empty_allowed, k, bLiteral);

        for (int idx = *iter; j < idx; ++j) {
            incr_j = true;
            n -= (m + 1);
            mpzIdx += temp;
            Counter->GetCount(temp, n, m, empty_allowed, k, false);
        }

        if (incr_j || (i + 1) >= (width - k)) {
            ++j;
            n -= m;
        }
    }
}

void rankPartsDistinctCapGmp(std::vector<int>::iterator iter,
                             int n, int m, int cap, int k,
                             double &dblIdx, mpz_class &mpzIdx) {

    const int width = m;
    mpzIdx = 0;

    n -= m;
    --cap;
    --m;

    mpz_class temp;
    const PartitionType ptype = PartitionType::DstctCapped;
    std::unique_ptr<CountClass> Counter = MakeCount(ptype);

    Counter->SetArrSize(ptype, n, m);
    Counter->InitializeMpz();

    std::vector<int> allowed(cap);
    std::iota(allowed.begin(), allowed.end(), 1);

    for (int i = 0, j = 0; i < (width - 1); ++i, n -= m, --m, ++j, ++iter) {
        Counter->GetCount(temp, n, m, allowed, k);

        for (int idx = *iter; j < idx; ++j) {
            n -= (m + 1);
            allowed.pop_back();
            mpzIdx += temp;
            Counter->GetCount(temp, n, m, allowed, k);
        }

        allowed.pop_back();
    }
}

void rankPartsDistinctCapMZGmp(std::vector<int>::iterator iter,
                               int n, int m, int cap, int k,
                               double &dblIdx, mpz_class &mpzIdx) {

    const int width = m;
    mpzIdx = 0;

    bool incr_j = false;
    --m;

    mpz_class temp;
    const PartitionType ptype = PartitionType::DstctCappedMZ;
    std::unique_ptr<CountClass> Counter = MakeCount(ptype);

    Counter->SetArrSize(ptype, n, m);
    Counter->InitializeMpz();

    std::vector<int> allowed(cap);
    std::iota(allowed.begin(), allowed.end(), 1);

    for (int i = 0, j = 0; i < (width - 1); ++i, --m, ++iter) {
        const bool bLiteral = !(incr_j || i >= (width - k));
        Counter->GetCount(temp, n, m, allowed, k, bLiteral);

        for (int idx = *iter; j < idx; ++j) {
            incr_j = true;
            n -= (m + 1);
            allowed.pop_back();
            mpzIdx += temp;
            Counter->GetCount(temp, n, m, allowed, k, false);
        }

        if (incr_j || (i + 1) >= (width - k)) {
            ++j;
            n -= m;
            allowed.pop_back();
        }
    }
}

rankPartsPtr GetRankPartsFunc(PartitionType ptype, bool IsGmp) {

    if (IsGmp) {
        switch (ptype) {
            case PartitionType::LengthOne: {
                return(rankPartsPtr(rankLengthOne));
            } case PartitionType::DstctCapped: {
                return(rankPartsPtr(rankPartsDistinctCapGmp));
            } case PartitionType::DstctCappedMZ: {
                return(rankPartsPtr(rankPartsDistinctCapMZGmp));
            } case PartitionType::DstctNoZero: {
                return(rankPartsPtr(rankPartsDistinctLenGmp));
            } case PartitionType::DstctOneZero: {
                return(rankPartsPtr(rankPartsDistinctOneZeroGmp));
            } case PartitionType::DstctMultiZero: {
                return(rankPartsPtr(rankPartsDistinctMultiZeroGmp));
            } case PartitionType::DstctStdAll: {
                return(rankPartsPtr(rankPartsDistinctMultiZeroGmp));
            } case PartitionType::RepCapped: {
                return(rankPartsPtr(rankPartsRepCapGmp));
            } case PartitionType::RepNoZero: {
                return(rankPartsPtr(rankPartsRepLenGmp));
            } case PartitionType::RepShort: {
                return(rankPartsPtr(rankPartsRepShortGmp));
            } case PartitionType::RepStdAll: {
                return(rankPartsPtr(rankPartsRepGmp));
            } case PartitionType::CompRepNoZero: {
                return(rankPartsPtr(rankCompsRepGmp));
            } case PartitionType::CompRepWeak: {
                return(rankPartsPtr(rankCompsRepGmp));
            } case PartitionType::CmpRpZroNotWk: {
                return(rankPartsPtr(rankCompsRepZeroGmp));
            } case PartitionType::CmpDstctNoZero: {
                return(rankPartsPtr(rankCompsDistinctGmp));
            } case PartitionType::CmpDstctCapped: {
                return(rankPartsPtr(rankCompsDistinctGmp));
            } case PartitionType::CmpDstctWeak: {
                return(rankPartsPtr(rankCompsDistinctGmp));
            } case PartitionType::CmpDstCapWeak: {
                return(rankPartsPtr(rankCompsDistinctGmp));
            } case PartitionType::CmpDstctZNotWk: {
                return(rankPartsPtr(rankCompsDistinctMZGmp));
            } case PartitionType::CmpDstCapMZNotWk: {
                return(rankPartsPtr(rankCompsDistinctMZGmp));
            } case PartitionType::CmpDstctMZWeak: {
                return(rankPartsPtr(rankCompsDistinctMZWeakGmp));
            } case PartitionType::CmpDstCapMZWeak: {
                return(rankPartsPtr(rankCompsDistinctMZWeakGmp));
            } default : {
                cpp11::stop("No algorithm available");
            }
        }
    } else {
        switch (ptype) {
            case PartitionType::LengthOne: {
                return(rankPartsPtr(rankLengthOne));
            } case PartitionType::DstctCapped: {
                return(rankPartsPtr(rankPartsDistinctCap));
            } case PartitionType::DstctCappedMZ: {
                return(rankPartsPtr(rankPartsDistinctCapMZ));
            } case PartitionType::DstctNoZero: {
                return(rankPartsPtr(rankPartsDistinctLen));
            } case PartitionType::DstctOneZero: {
                return(rankPartsPtr(rankPartsDistinctOneZero));
            } case PartitionType::DstctMultiZero: {
                return(rankPartsPtr(rankPartsDistinctMultiZero));
            } case PartitionType::DstctStdAll: {
                return(rankPartsPtr(rankPartsDistinctMultiZero));
            } case PartitionType::RepCapped: {
                return(rankPartsPtr(rankPartsRepCap));
            } case PartitionType::RepNoZero: {
                return(rankPartsPtr(rankPartsRepLen));
            } case PartitionType::RepShort: {
                return(rankPartsPtr(rankPartsRepShort));
            } case PartitionType::RepStdAll: {
                return(rankPartsPtr(rankPartsRep));
            } case PartitionType::CompRepNoZero: {
                return(rankPartsPtr(rankCompsRep));
            } case PartitionType::CompRepWeak: {
                return(rankPartsPtr(rankCompsRep));
            } case PartitionType::CmpRpZroNotWk: {
                return(rankPartsPtr(rankCompsRepZero));
            } case PartitionType::CmpDstctNoZero: {
                return(rankPartsPtr(rankCompsDistinct));
            } case PartitionType::CmpDstctCapped: {
                return(rankPartsPtr(rankCompsDistinct));
            } case PartitionType::CmpDstctWeak: {
                return(rankPartsPtr(rankCompsDistinct));
            } case PartitionType::CmpDstCapWeak: {
                return(rankPartsPtr(rankCompsDistinct));
            } case PartitionType::CmpDstctZNotWk: {
                return(rankPartsPtr(rankCompsDistinctMZ));
            } case PartitionType::CmpDstCapMZNotWk: {
                return(rankPartsPtr(rankCompsDistinctMZ));
            } case PartitionType::CmpDstctMZWeak: {
                return(rankPartsPtr(rankCompsDistinctMZWeak));
            } case PartitionType::CmpDstCapMZWeak: {
                return(rankPartsPtr(rankCompsDistinctMZWeak));
            } default : {
                cpp11::stop("No algorithm available");
            }
        }
    }
}
