#include "Permutations/BigPermuteCount.h"
#include "Permutations/PermuteCount.h"
#include <algorithm> // std::find
#include <numeric>
#include <cmath>

using nthPermPtr = std::vector<int> (*const)(int n, int m, double dblIdx,
                                     const mpz_class &mpzIdx,
                                     const std::vector<int> &Reps);

std::vector<int> nthPerm(int n, int m, double dblIdx, const mpz_class &mpzIdx,
                         const std::vector<int> &Reps) {

    double index1 = dblIdx;
    std::vector<int> res(m);

    double temp = NumPermsNoRep(n, m);
    std::vector<int> indexVec(n);
    std::iota(indexVec.begin(), indexVec.end(), 0);

    for (int k = 0, n1 = n; k < m; ++k, --n1) {
        temp /= n1;
        int j = static_cast<int>(index1 / temp);
        res[k] = indexVec[j];
        index1 -= (temp * j);
        indexVec.erase(indexVec.begin() + j);
    }

    return res;
}

std::vector<int> nthPermRep(int n, int m, double dblIdx, const mpz_class &mpzIdx,
                            const std::vector<int> &Reps) {

    double index1 = dblIdx;
    std::vector<int> res(m);

    double temp = std::pow(static_cast<double>(n),
                           static_cast<double>(m));

    for (int k = 0; k < m; ++k) {
        temp /= n;
        int j = static_cast<int>(index1 / temp);
        res[k] = j;
        index1 -= (temp * j);
    }

    return res;
}

std::vector<int> nthPermMult(int n, int m, double dblIdx, const mpz_class &mpzIdx,
                             const std::vector<int> &Reps) {

    double index1 = dblIdx + 1;
    double index2 = index1;

    std::vector<int> res(m);
    std::vector<int> Counts;
    std::vector<int> TempReps = Reps;

    for (int k = 0, r1 = m - 1; k < m; ++k, --r1) {

        int j = 0;

        while (TempReps[j] == 0) {
            ++j;
        }

        --TempReps[j];
        Counts = nonZeroVec(TempReps);
        double test = MultisetPermRowNum(Counts.size(), r1, Counts);
        double temp = test;

        for (; test < index1; test += temp) {
            index2 -= temp;
            ++TempReps[j];
            ++j;

            while (TempReps[j] == 0) {
                ++j;
            }

            --TempReps[j];
            Counts = nonZeroVec(TempReps);
            temp = MultisetPermRowNum(Counts.size(), r1, Counts);
        }

        res[k] = j;
        index1 = index2;
    }

    return res;
}

std::vector<int> nthPermGmp(int n, int m, double dblIdx,
                            const mpz_class &mpzIdx,
                            const std::vector<int> &Reps) {

    mpz_class temp;
    mpz_class temp2;
    mpz_class index1(mpzIdx);
    NumPermsNoRepGmp(temp, n, m);

    std::vector<int> res(m);
    std::vector<int> indexVec(n);
    std::iota(indexVec.begin(), indexVec.end(), 0);

    for (int k = 0, n1 = n; k < m; ++k, --n1) {
        mpz_divexact_ui(temp.get_mpz_t(), temp.get_mpz_t(), n1);
        temp2 = index1 / temp;
        int j = temp2.get_si();
        res[k] = indexVec[j];
        index1 -= (temp * j);
        indexVec.erase(indexVec.begin() + j);
    }

    return res;
}

std::vector<int> nthPermRepGmp(int n, int m, double dblIdx,
                               const mpz_class &mpzIdx,
                               const std::vector<int> &Reps) {

    mpz_class temp;
    mpz_class temp2;
    mpz_class index1(mpzIdx);

    std::vector<int> res(m);
    mpz_ui_pow_ui(temp.get_mpz_t(), n, m);

    for (int k = 0; k < m; ++k) {
        mpz_divexact_ui(temp.get_mpz_t(), temp.get_mpz_t(), n);
        temp2 = index1 / temp;
        int j = temp2.get_si();
        res[k] = j;
        index1 -= (temp * j);
    }

    return res;
}

std::vector<int> nthPermMultGmp(int n, int m, double dblIdx,
                                const mpz_class &mpzIdx,
                                const std::vector<int> &Reps) {

    mpz_class temp;
    mpz_class index1(mpzIdx);
    ++index1;

    std::vector<int> res(m);
    std::vector<int> Counts;
    std::vector<int> TempReps = Reps;

    mpz_class test;
    mpz_class index2(index1);

    for (int k = 0, r1 = m - 1; k < m; ++k, --r1) {

        int j = 0;

        while (TempReps[j] == 0) {
            ++j;
        }

        --TempReps[j];
        Counts = nonZeroVec(TempReps);
        MultisetPermRowNumGmp(temp, static_cast<int>(Counts.size()),
                              r1, Counts);
        test = temp;

        while (cmp(test, index1) < 0) {
            index2 -= temp;
            ++TempReps[j];
            ++j;

            while (TempReps[j] == 0) {
                ++j;
            }

            --TempReps[j];
            Counts = nonZeroVec(TempReps);
            MultisetPermRowNumGmp(temp, static_cast<int>(Counts.size()),
                                  r1, Counts);
            test += temp;
        }

        res[k] = j;
        index1 = index2;
    }

    return res;
}

nthPermPtr GetNthPermFunc(bool IsMult, bool IsRep, bool IsGmp) {

    if (IsGmp) {
        if (IsMult) {
            return(nthPermPtr(nthPermMultGmp));
        } else if (IsRep) {
            return(nthPermPtr(nthPermRepGmp));
        } else {
            return(nthPermPtr(nthPermGmp));
        }
    } else {
        if (IsMult) {
            return(nthPermPtr(nthPermMult));
        } else if (IsRep) {
            return(nthPermPtr(nthPermRep));
        } else {
            return(nthPermPtr(nthPerm));
        }
    }
}

void TopOffPerm(std::vector<int> &z, const std::vector<int> &myReps,
                int n, int m, bool IsRep, bool IsMult) {

    if (IsMult) {
        std::vector<int> f(n, 0);

        for (int i = 0; i < m; ++i) {
            ++f[z[i]];
        }

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < (myReps[i] - f[i]); ++j) {
                z.push_back(i);
            }
        }
    } else if (!IsRep && m < n) {
        for (int i = 0; i < n; ++i) {
            if (std::find(z.begin(), z.end(), i) == z.end()) {
                z.push_back(i);
            }
        }
    }
}
