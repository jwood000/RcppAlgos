#include "Partitions/RankPartition.h"
#include "Partitions/PartitionsCountDistinct.h"
#include "Partitions/PartitionsCountRep.h"
#include "Partitions/PartitionsCount.h"
#include <memory>
#include <numeric>
#include <vector>

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

void rankPartsRepShort(std::vector<int>::iterator iter, int n, int m,
                       int cap, int k, double &dblIdx, mpz_class &mpzIdx) {

    rankPartsRepLen(iter, n, m, cap, k, dblIdx, mpzIdx);
}

void rankPartsRepShortGmp(std::vector<int>::iterator iter, int n, int m,
                          int cap, int k, double &dblIdx, mpz_class &mpzIdx) {

    rankPartsRepLenGmp(iter, n, m, cap, k, dblIdx, mpzIdx);
}

void rankPartsRep(std::vector<int>::iterator iter, int n, int m,
                  int cap, int k, double &dblIdx, mpz_class &mpzIdx) {

    rankPartsRepLen(iter, n * 2, m, cap, k, dblIdx, mpzIdx);
}

void rankPartsRepGmp(std::vector<int>::iterator iter, int n, int m,
                     int cap, int k, double &dblIdx, mpz_class &mpzIdx) {

    rankPartsRepLenGmp(iter, n * 2, m, cap, k, dblIdx, mpzIdx);
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

void rankPartsDistinctLenGmp(std::vector<int>::iterator iter, int n, int m,
                            int cap, int k, double &dblIdx, mpz_class &mpzIdx) {

    const int width = m;
    mpzIdx = 0;

    n -= m;
    --m;

    mpz_class temp;
    const PartitionType ptype = PartitionType::DistinctNoZero;
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

void rankPartsDistinctOneZero(
    std::vector<int>::iterator iter, int n, int m,
    int cap, int k, double &dblIdx, mpz_class &mpzIdx
) {

    rankPartsDistinctLen(iter, n, m, cap, k, dblIdx, mpzIdx);
}

void rankPartsDistinctOneZeroGmp(std::vector<int>::iterator iter,
                                 int n, int m, int cap, int k,
                                 double &dblIdx, mpz_class &mpzIdx) {

    rankPartsDistinctLenGmp(iter, n, m, cap, k, dblIdx, mpzIdx);
}

void rankPartsDistinctMZ(
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
                      CountPartsDistinctMZ(n, m, empty_allowed, k);

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

void rankPartsDistinctMZGmp(std::vector<int>::iterator iter,
                                   int n, int m, int cap, int k,
                                   double &dblIdx, mpz_class &mpzIdx) {

    const int width = m;
    mpzIdx = 0;

    bool incr_j = false;
    --m;

    mpz_class temp;
    const PartitionType ptype = PartitionType::DistinctMZ;
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

void rankPartsDistinctCapGmp(std::vector<int>::iterator iter,
                             int n, int m, int cap, int k,
                             double &dblIdx, mpz_class &mpzIdx) {

    const int width = m;
    mpzIdx = 0;

    n -= m;
    --cap;
    --m;

    mpz_class temp;
    const PartitionType ptype = PartitionType::DistinctCapped;
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

void rankPartsDistinctCapMZGmp(std::vector<int>::iterator iter,
                               int n, int m, int cap, int k,
                               double &dblIdx, mpz_class &mpzIdx) {

    const int width = m;
    mpzIdx = 0;

    bool incr_j = false;
    --m;

    mpz_class temp;
    const PartitionType ptype = PartitionType::DistinctCappedMZ;
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
