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
using rankPartsPtr = void (*const)(std::vector<int>::iterator iter,
                           int n, int m, int cap, int k,
                           double &dblIdx, mpz_class &mpzIdx);

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
    const int max_val = n - (width * (width - 1)) / 2;
    dblIdx = 0;

    std::vector<char> mask(n + 1, 0);
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
            while (mask[j + 1]) {
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

        while (mask[j + 1]) {
            ++j;
        }

        cur_val = j + 1;
        partial_sum += (j + 1);
        UpdateAllowed(mask, allowed, i + 1, j + 1, width,
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
    const int max_val = n - (width * (width - 1)) / 2;
    mpzIdx = 0;

    std::vector<char> mask(n + 1, 0);
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
            while (mask[j + 1]) {
                ++j;
            }

            partial_sum += (j + 1 - cur_val);
            UpdateAllowed(mask, allowed, i, j + 1, width,
                          n, cur_val, partial_sum);

            mpzIdx += temp;
            temp = CountCompDistLenRstrctd(n - partial_sum, m, allowed);

            ++j;
            keepGoing = j <= idx;
        }

        j = 0;

        while (mask[j + 1]) {
            ++j;
        }

        cur_val = j + 1;
        partial_sum += (j + 1);
        UpdateAllowed(mask, allowed, i + 1, j + 1, width,
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
            } default : {
                cpp11::stop("No algorithm available");
            }
        }
    }
}
