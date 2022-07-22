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
                           double &dblIdx, mpz_t mpzIdx);

//*********************** Compositions Funcitons **************************//

void rankCompsRep(std::vector<int>::iterator iter, int n, int m,
                  int cap, int k, double &dblIdx, mpz_t mpzIdx) {

    const int width = m;
    dblIdx = 0;

    --n;
    --m;

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m, j = 0, ++iter) {
        double temp = CountCompsRepLen(n, m, cap, k);

        for (int idx = *iter; j < idx; ++j) {
            --n;
            dblIdx += temp;
            temp = CountCompsRepLen(n, m, cap, k);
        }
    }
}

void rankCompsRepZero(std::vector<int>::iterator iter, int n, int m,
                      int cap, int k, double &dblIdx, mpz_t mpzIdx) {

    const int width = m;
    dblIdx = 0;

    bool incr_j = false;
    --m;

    for (int i = 0, j = 0; i < (width - 1); ++i, --m, j = incr_j, ++iter) {
        double temp = incr_j ? CountCompsRepLen(n, m, cap, k) :
                               CountCompsRepZero(n, m, cap, k);

        for (int idx = *iter; j < idx; ++j) {
            incr_j = true;
            --n;
            dblIdx += temp;
            temp = CountCompsRepLen(n, m, cap, k);
        }

        if (incr_j) --n;
    }
}

//************************* Paritions Funcitons ***************************//

void rankPartsRepLen(std::vector<int>::iterator iter, int n, int m,
                     int cap, int k, double &dblIdx, mpz_t mpzIdx) {

    const int width = m;
    dblIdx = 0;

    --n;
    --m;

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m, ++iter) {
        double temp = CountPartsRepLen(n, m, cap, k);

        for (int idx = *iter; j < idx; ++j) {
            n -= (m + 1);
            dblIdx += temp;
            temp = CountPartsRepLen(n, m, cap, k);
        }
    }
}

void rankPartsRepShort(std::vector<int>::iterator iter, int n, int m,
                       int cap, int k, double &dblIdx, mpz_t mpzIdx) {

    rankPartsRepLen(iter, n, m, cap, k, dblIdx, mpzIdx);
}

void rankPartsRep(std::vector<int>::iterator iter, int n, int m,
                  int cap, int k, double &dblIdx, mpz_t mpzIdx) {

    rankPartsRepLen(iter, n * 2, m, cap, k, dblIdx, mpzIdx);
}

void rankPartsRepCap(std::vector<int>::iterator iter, int n, int m,
                     int cap, int k, double &dblIdx, mpz_t mpzIdx) {

    const int width = m;
    dblIdx = 0;

    --n;
    --m;

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m, ++iter) {
        double temp = CountPartsRepLenCap(n, m, cap, k);

        for (int idx = *iter; j < idx; ++j) {
            n -= (m + 1);
            --cap;
            dblIdx += temp;
            temp = CountPartsRepLenCap(n, m, cap, k);
        }
    }
}

void rankPartsDistinctLen(std::vector<int>::iterator iter, int n, int m,
                          int cap, int k, double &dblIdx, mpz_t mpzIdx) {

    const int width = m;
    dblIdx = 0;

    n -= m;
    --m;

    for (int i = 0, j = 0; i < (width - 1); ++i, n -= m, --m, ++j, ++iter) {
        double temp = CountPartsDistinctLen(n, m, cap, k);

        for (int idx = *iter; j < idx; ++j) {
            n -= (m + 1);
            dblIdx += temp;
            temp = CountPartsDistinctLen(n, m, cap, k);
        }
    }
}

void rankPartsDistinctOneZero(std::vector<int>::iterator iter, int n, int m,
                              int cap, int k, double &dblIdx, mpz_t mpzIdx) {

    rankPartsDistinctLen(iter, n, m, cap, k, dblIdx, mpzIdx);
}

void rankPartsDistinctMultiZero(std::vector<int>::iterator iter, int n, int m,
                                int cap, int k, double &dblIdx, mpz_t mpzIdx) {

    const int width = m;
    dblIdx = 0;

    bool incr_j = false;
    --m;

    for (int i = 0, j = 0; i < (width - 1); ++i, --m, ++iter) {
        double temp = (incr_j || i >= (width - k)) ?
                      CountPartsDistinctLen(n, m, cap, k) :
                      CountPartsDistinctMultiZero(n, m, cap, k);

        for (int idx = *iter; j < idx; ++j) {
            incr_j = true;
            n -= (m + 1);
            dblIdx += temp;
            temp = CountPartsDistinctLen(n, m, cap, k);
        }

        if (incr_j || (i + 1) >= (width - k)) {
            ++j;
            n -= m;
        }
    }
}

void rankPartsDistinctCap(std::vector<int>::iterator iter, int n, int m,
                          int cap, int k, double &dblIdx, mpz_t mpzIdx) {

    const int width = m;
    dblIdx = 0;

    n -= m;
    --cap;
    --m;

    for (int i = 0, j = 0; i < (width - 1);
         ++i, n -= m, --m, ++j, --cap, ++iter) {
        double temp = CountPartsDistinctLenCap(n, m, cap, k);

        for (int idx = *iter; j < idx; ++j) {
            n -= (m + 1);
            --cap;
            dblIdx += temp;
            temp = CountPartsDistinctLenCap(n, m, cap, k);
        }
    }
}

void rankPartsDistinctCapMZ(std::vector<int>::iterator iter, int n, int m,
                            int cap, int k, double &dblIdx, mpz_t mpzIdx) {

    const int width = m;
    dblIdx = 0;

    bool incr_j = false;
    --m;

    for (int i = 0, j = 0; i < (width - 1); ++i, --m, ++iter) {
        double temp = (incr_j || i >= (width - k)) ?
                      CountPartsDistinctLenCap(n, m, cap, k) :
                      CountPartsDistinctCapMZ(n, m, cap, k);

        for (int idx = *iter; j < idx; ++j) {
            incr_j = true;
            n -= (m + 1);
            --cap;
            dblIdx += temp;
            temp = CountPartsDistinctLenCap(n, m, cap, k);
        }

        if (incr_j || (i + 1) >= (width - k)) {
            ++j;
            n -= m;
            --cap;
        }
    }
}

//*********************** Starting Gmp Funcitons **************************//

void rankCompsRepGmp(std::vector<int>::iterator iter, int n, int m,
                     int cap, int k, double &dblIdx, mpz_t mpzIdx) {

    const int width = m;
    mpz_set_ui(mpzIdx, 0u);

    --n;
    --m;

    mpz_t temp;
    mpz_init(temp);

    const PartitionType ptype = PartitionType::RepNoZero;
    std::unique_ptr<CountClass> myClass = MakeCount(ptype, true);

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m, j = 0, ++iter) {
        myClass->GetCount(temp, n, m, cap, k);

        for (int idx = *iter; j < idx; ++j) {
            --n;
            mpz_add(mpzIdx, mpzIdx, temp);
            myClass->GetCount(temp, n, m, cap, k);
        }
    }

    mpz_clear(temp);
}

void rankCompsRepZeroGmp(std::vector<int>::iterator iter, int n, int m,
                         int cap, int k, double &dblIdx, mpz_t mpzIdx) {

    const int width = m;
    mpz_set_ui(mpzIdx, 0u);

    bool incr_j = false;
    --m;

    mpz_t temp;
    mpz_init(temp);

    const PartitionType ptype = PartitionType::RepShort;
    std::unique_ptr<CountClass> myClass = MakeCount(ptype, true);

    for (int i = 0, j = 0; i < (width - 1); ++i, --m, j = incr_j, ++iter) {
        myClass->GetCount(temp, n, m, cap, k, !incr_j);

        for (int idx = *iter; j < idx; ++j) {
            incr_j = true;
            --n;
            mpz_add(mpzIdx, mpzIdx, temp);
            myClass->GetCount(temp, n, m, cap, k, false);
        }

        mpz_set_ui(temp, 0);
        if (incr_j) --n;
    }

    mpz_clear(temp);
}

void rankPartsRepLenGmp(std::vector<int>::iterator iter, int n, int m,
                        int cap, int k, double &dblIdx, mpz_t mpzIdx) {

    const int width = m;
    mpz_set_ui(mpzIdx, 0u);

    --n;
    --m;

    mpz_t temp;
    mpz_init(temp);

    const PartitionType ptype = PartitionType::RepShort;
    std::unique_ptr<CountClass> myClass = MakeCount(ptype);

    myClass->SetArrSize(ptype, n, m, cap);
    myClass->InitializeMpz();

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m, ++iter) {
        myClass->GetCount(temp, n, m, cap, k);

        for (int idx = *iter; j < idx; ++j) {
            n -= (m + 1);
            mpz_add(mpzIdx, mpzIdx, temp);
            myClass->GetCount(temp, n, m, cap, k);
        }
    }

    mpz_clear(temp);
    myClass->ClearMpz();
}

void rankPartsRepShortGmp(std::vector<int>::iterator iter, int n, int m,
                          int cap, int k, double &dblIdx, mpz_t mpzIdx) {

    rankPartsRepLenGmp(iter, n, m, cap, k, dblIdx, mpzIdx);
}

void rankPartsRepGmp(std::vector<int>::iterator iter, int n, int m,
                     int cap, int k, double &dblIdx, mpz_t mpzIdx) {

    rankPartsRepLenGmp(iter, n * 2, m, cap, k, dblIdx, mpzIdx);
}

void rankPartsRepCapGmp(std::vector<int>::iterator iter, int n, int m,
                        int cap, int k, double &dblIdx, mpz_t mpzIdx) {

    const int width = m;
    mpz_set_ui(mpzIdx, 0u);

    --n;
    --m;

    mpz_t temp;
    mpz_init(temp);

    const PartitionType ptype = PartitionType::RepCapped;
    std::unique_ptr<CountClass> myClass = MakeCount(ptype);

    myClass->SetArrSize(ptype, n, m, cap);
    myClass->InitializeMpz();

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m, ++iter) {
        myClass->GetCount(temp, n, m, cap, k);

        for (int idx = *iter; j < idx; ++j) {
            n -= (m + 1);
            --cap;
            mpz_add(mpzIdx, mpzIdx, temp);
            myClass->GetCount(temp, n, m, cap, k);
        }
    }

    mpz_clear(temp);
    myClass->ClearMpz();
}

void rankPartsDistinctLenGmp(std::vector<int>::iterator iter, int n, int m,
                             int cap, int k, double &dblIdx, mpz_t mpzIdx) {

    const int width = m;
    mpz_set_ui(mpzIdx, 0u);

    n -= m;
    --m;

    mpz_t temp;
    mpz_init(temp);

    const PartitionType ptype = PartitionType::DstctNoZero;
    std::unique_ptr<CountClass> myClass = MakeCount(ptype);

    myClass->SetArrSize(ptype, n, m, cap);
    myClass->InitializeMpz();

    for (int i = 0, j = 0; i < (width - 1); ++i, n -= m, --m, ++j, ++iter) {
        myClass->GetCount(temp, n, m, cap, k);

        for (int idx = *iter; j < idx; ++j) {
            n -= (m + 1);
            mpz_add(mpzIdx, mpzIdx, temp);
            myClass->GetCount(temp, n, m, cap, k);
        }
    }

    mpz_clear(temp);
    myClass->ClearMpz();
}

void rankPartsDistinctOneZeroGmp(std::vector<int>::iterator iter,
                                 int n, int m, int cap, int k,
                                 double &dblIdx, mpz_t mpzIdx) {

    rankPartsDistinctLenGmp(iter, n, m, cap, k, dblIdx, mpzIdx);
}

void rankPartsDistinctMultiZeroGmp(std::vector<int>::iterator iter,
                                   int n, int m, int cap, int k,
                                   double &dblIdx, mpz_t mpzIdx) {

    const int width = m;
    mpz_set_ui(mpzIdx, 0u);

    bool incr_j = false;
    --m;

    mpz_t temp;
    mpz_init(temp);

    const PartitionType ptype = PartitionType::DstctMultiZero;
    std::unique_ptr<CountClass> myClass = MakeCount(ptype);

    myClass->SetArrSize(ptype, n, m, cap);
    myClass->InitializeMpz();

    for (int i = 0, j = 0; i < (width - 1); ++i, --m, ++iter) {
        const bool bLiteral = !(incr_j || i >= (width - k));
        myClass->GetCount(temp, n, m, cap, k, bLiteral);

        for (int idx = *iter; j < idx; ++j) {
            incr_j = true;
            n -= (m + 1);
            mpz_add(mpzIdx, mpzIdx, temp);
            myClass->GetCount(temp, n, m, cap, k, false);
        }

        if (incr_j || (i + 1) >= (width - k)) {
            ++j;
            n -= m;
        }
    }

    mpz_clear(temp);
    myClass->ClearMpz();
}

void rankPartsDistinctCapGmp(std::vector<int>::iterator iter,
                             int n, int m, int cap, int k,
                             double &dblIdx, mpz_t mpzIdx) {

    const int width = m;
    mpz_set_ui(mpzIdx, 0u);

    n -= m;
    --cap;
    --m;

    mpz_t temp;
    mpz_init(temp);

    const PartitionType ptype = PartitionType::DstctCapped;
    std::unique_ptr<CountClass> myClass = MakeCount(ptype);

    myClass->SetArrSize(ptype, n, m, cap);
    myClass->InitializeMpz();

    for (int i = 0, j = 0; i < (width - 1);
         ++i, n -= m, --m, ++j, --cap, ++iter) {
        myClass->GetCount(temp, n, m, cap, k);

        for (int idx = *iter; j < idx; ++j) {
            n -= (m + 1);
            --cap;
            mpz_add(mpzIdx, mpzIdx, temp);
            myClass->GetCount(temp, n, m, cap, k);
        }
    }

    mpz_clear(temp);
    myClass->ClearMpz();
}

void rankPartsDistinctCapMZGmp(std::vector<int>::iterator iter,
                               int n, int m, int cap, int k,
                               double &dblIdx, mpz_t mpzIdx) {

    const int width = m;
    mpz_set_ui(mpzIdx, 0u);

    bool incr_j = false;
    --m;

    mpz_t temp;
    mpz_init(temp);

    const PartitionType ptype = PartitionType::DstctCappedMZ;
    std::unique_ptr<CountClass> myClass = MakeCount(ptype);

    myClass->SetArrSize(ptype, n, m, cap);
    myClass->InitializeMpz();

    for (int i = 0, j = 0; i < (width - 1); ++i, --m, ++iter) {
        const bool bLiteral = !(incr_j || i >= (width - k));
        myClass->GetCount(temp, n, m, cap, k, bLiteral);

        for (int idx = *iter; j < idx; ++j) {
            incr_j = true;
            n -= (m + 1);
            --cap;
            mpz_add(mpzIdx, mpzIdx, temp);
            myClass->GetCount(temp, n, m, cap, k, false);
        }

        if (incr_j || (i + 1) >= (width - k)) {
            ++j;
            n -= m;
            --cap;
        }
    }

    mpz_clear(temp);
    myClass->ClearMpz();
}

rankPartsPtr GetRankPartsFunc(PartitionType ptype, bool IsGmp, bool IsComp) {

    if (IsComp && IsGmp) {
        switch (ptype) {
            case PartitionType::RepNoZero : {
                return(rankPartsPtr(rankCompsRepGmp));
            } case PartitionType::RepShort : {
                return(rankPartsPtr(rankCompsRepZeroGmp));
            } case PartitionType::RepStdAll : {
                return(rankPartsPtr(rankCompsRepZeroGmp));
            }default : {
                cpp11::stop("No algorithm available");
            }
        }
    } else if (IsComp) {
        switch (ptype) {
            case PartitionType::RepNoZero : {
                return(rankPartsPtr(rankCompsRep));
            } case PartitionType::RepShort : {
                return(rankPartsPtr(rankCompsRepZero));
            } case PartitionType::RepStdAll : {
                return(rankPartsPtr(rankCompsRepZero));
            } default : {
                cpp11::stop("No algorithm available");
            }
        }
    } else if (IsGmp) {
        switch (ptype) {
            case PartitionType::DstctCapped: {
                return(rankPartsPtr(rankPartsDistinctCapGmp));
            } case PartitionType::DstctCappedMZ: {
                return(rankPartsPtr(rankPartsDistinctCapMZGmp));
            } case PartitionType::DstctNoZero : {
                return(rankPartsPtr(rankPartsDistinctLenGmp));
            } case PartitionType::DstctOneZero: {
                return(rankPartsPtr(rankPartsDistinctOneZeroGmp));
            } case PartitionType::DstctMultiZero : {
                return(rankPartsPtr(rankPartsDistinctMultiZeroGmp));
            } case PartitionType::DstctStdAll: {
                return(rankPartsPtr(rankPartsDistinctMultiZeroGmp));
            } case PartitionType::RepCapped : {
                return(rankPartsPtr(rankPartsRepCapGmp));
            } case PartitionType::RepNoZero: {
                return(rankPartsPtr(rankPartsRepLenGmp));
            } case PartitionType::RepShort : {
                return(rankPartsPtr(rankPartsRepShortGmp));
            } case PartitionType::RepStdAll : {
                return(rankPartsPtr(rankPartsRepGmp));
            } default : {
                cpp11::stop("No algorithm available");
            }
        }
    } else {
        switch (ptype) {
            case PartitionType::DstctCapped: {
                return(rankPartsPtr(rankPartsDistinctCap));
            } case PartitionType::DstctCappedMZ: {
                return(rankPartsPtr(rankPartsDistinctCapMZ));
            } case PartitionType::DstctNoZero : {
                return(rankPartsPtr(rankPartsDistinctLen));
            } case PartitionType::DstctOneZero: {
                return(rankPartsPtr(rankPartsDistinctOneZero));
            } case PartitionType::DstctMultiZero : {
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
            } default : {
                cpp11::stop("No algorithm available");
            }
        }
    }
}
