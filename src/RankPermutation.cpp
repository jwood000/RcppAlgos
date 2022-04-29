#include "Permutations/BigPermuteCount.h"
#include "Permutations/PermuteCount.h"
#include <algorithm> // std::find
#include <numeric>
#include <cmath>

using It = std::vector<int>::iterator;

using rankPermPtr = void (*const)(It iter, int n, int r, double &dblIdx,
                          mpz_t mpzIdx, const std::vector<int> &Reps);

int which(const std::vector<int> &idx, int j) {
    auto it = std::find(idx.cbegin(), idx.cend(), j);
    return std::distance(idx.cbegin(), it);
}

void rankPerm(It iter, int n, int r, double &dblIdx,
              mpz_t mpzIdx, const std::vector<int> &Reps) {

    dblIdx = 0;
    double temp = NumPermsNoRep(n, r);

    std::vector<int> indexVec(n);
    std::iota(indexVec.begin(), indexVec.end(), 0);

    for (int k = 0, n1 = n; k < r; ++k, --n1, ++iter) {
        temp /= n1;
        int j = which(indexVec, *iter);
        dblIdx += (temp * j);
        indexVec.erase(indexVec.begin() + j);
    }
}

void rankPermRep(It iter, int n, int r, double &dblIdx,
                 mpz_t mpzIdx, const std::vector<int> &Repss) {

    dblIdx = 0;
    double temp = std::pow(static_cast<double>(n),
                           static_cast<double>(r));

    for (int k = 0; k < r; ++k, ++iter) {
        temp /= n;
        dblIdx += (temp * (*iter));
    }
}

void rankPermMult(It iter, int n, int r, double &dblIdx,
                  mpz_t mpzIdx, const std::vector<int> &Reps) {

    dblIdx = 0;
    std::vector<int> Counts;
    std::vector<int> TempReps = Reps;

    for (int k = 0, r1 = r - 1; k < r; ++k, --r1, ++iter) {

        int j = 0;

        while (TempReps[j] == 0) {
            ++j;
        }

        --TempReps[j];
        Counts = nonZeroVec(TempReps);
        double test = MultisetPermRowNum(Counts.size(), r1, Counts);
        double temp = test;

        for (int idx = *iter; j < idx; test += temp) {
            dblIdx += temp;
            ++TempReps[j];
            ++j;

            while (TempReps[j] == 0) {
                ++j;
            }

            --TempReps[j];
            Counts = nonZeroVec(TempReps);
            temp = MultisetPermRowNum(Counts.size(), r1, Counts);
        }
    }
}

void rankPermGmp(It iter, int n, int r, double &dblIdx,
                 mpz_t mpzIdx, const std::vector<int> &Reps) {

    mpz_t temp;
    mpz_init(temp);

    mpz_set_ui(mpzIdx, 0u);
    NumPermsNoRepGmp(temp, n, r);

    std::vector<int> indexVec(n);
    std::iota(indexVec.begin(), indexVec.end(), 0);

    for (int k = 0, n1 = n; k < r; ++k, --n1, ++iter) {
        mpz_divexact_ui(temp, temp, n1);
        int j = which(indexVec, *iter);
        mpz_addmul_ui(mpzIdx, temp, j);
        indexVec.erase(indexVec.begin() + j);
    }

    mpz_clear(temp);
}

void rankPermRepGmp(It iter, int n, int r, double &dblIdx,
                    mpz_t mpzIdx, const std::vector<int> &Reps) {

    mpz_t temp;
    mpz_init(temp);
    mpz_set_ui(mpzIdx, 0u);
    mpz_ui_pow_ui(temp, n, r);

    for (int k = 0; k < r; ++k, ++iter) {
        mpz_divexact_ui(temp, temp, n);
        mpz_addmul_ui(mpzIdx, temp, *iter);
    }

    mpz_clear(temp);
}

void rankPermMultGmp(It iter, int n, int r, double &dblIdx,
                     mpz_t mpzIdx, const std::vector<int> &Reps) {

    mpz_t temp;
    mpz_init(temp);
    mpz_set_ui(mpzIdx, 0u);

    std::vector<int> Counts;
    std::vector<int> TempReps = Reps;

    mpz_t test;
    mpz_init(test);

    for (int k = 0, r1 = r - 1; k < r; ++k, --r1, ++iter) {

        int j = 0;

        while (TempReps[j] == 0) {
            ++j;
        }

        --TempReps[j];
        Counts = nonZeroVec(TempReps);
        MultisetPermRowNumGmp(temp, static_cast<int>(Counts.size()), r1, Counts);
        mpz_set(test, temp);

        for (int idx = *iter; j < idx;) {
            mpz_add(mpzIdx, mpzIdx, temp);
            ++TempReps[j];
            ++j;

            while (TempReps[j] == 0) {
                ++j;
            }

            --TempReps[j];
            Counts = nonZeroVec(TempReps);
            MultisetPermRowNumGmp(temp, static_cast<int>(Counts.size()), r1, Counts);
            mpz_add(test, test, temp);
        }
    }

    mpz_clear(test);
    mpz_clear(temp);
}

rankPermPtr GetRankPermFunc(bool IsMult, bool IsRep, bool IsGmp) {

    if (IsGmp) {
        if (IsMult) {
            return(rankPermPtr(rankPermMultGmp));
        } else if (IsRep) {
            return(rankPermPtr(rankPermRepGmp));
        } else {
            return(rankPermPtr(rankPermGmp));
        }
    } else {
        if (IsMult) {
            return(rankPermPtr(rankPermMult));
        } else if (IsRep) {
            return(rankPermPtr(rankPermRep));
        } else {
            return(rankPermPtr(rankPerm));
        }
    }
}
