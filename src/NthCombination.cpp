#include "Combinations/BigComboCount.h"
#include "Combinations/ComboCount.h"

using nthCombPtr = std::vector<int> (*const)(int n, int m, double dblIdx,
                                     const mpz_class &mpzIdx,
                                     const std::vector<int> &Reps);

std::vector<int> nthComb(int n, int m, double dblIdx,
                         const mpz_class &mpzIdx,
                         const std::vector<int> &Reps) {

    std::vector<int> res(m);
    double temp = nChooseK(n - 1, m - 1);

    for (int k = 0, j = 0, n1 = n - 1,
         r1 = m - 1; k < m; ++k, --n1, --r1, ++j) {

        for (int s = n1 - r1; temp <= dblIdx; --n1, ++j, --s) {
            dblIdx -= temp;
            temp *= s;
            temp /= n1;
        }

        temp *= r1;
        temp /= n1;
        res[k] = j;
    }

    return res;
}

std::vector<int> nthCombRep(int n, int m, double dblIdx,
                            const mpz_class &mpzIdx,
                            const std::vector<int> &Reps) {

    std::vector<int> res(m);
    double temp = NumCombsWithRep(n, m - 1);

    for (int k = 0, j = 0, n1 = n, r1 = m - 1; k < m; ++k, --r1) {
        for (; temp <= dblIdx; --n1, ++j) {
            dblIdx -= temp;
            temp *= (n1 - 1);
            temp /= (n1 + r1 - 1);
        }

        temp *= r1;
        temp /= (n1 + r1 - 1);
        res[k] = j;
    }

    return res;
}

std::vector<int> nthCombMult(int n, int m, double dblIdx,
                             const mpz_class &mpzIdx,
                             const std::vector<int> &Reps) {

    std::vector<int> res(m);
    std::vector<int> Counts = Reps;
    std::vector<int> TempReps = Reps;

    for (int k = 0, j = 0, n1 = n, r1 = m - 1; k < m; ++k, --r1) {

        --Counts.front();
        if (Counts.front() == 0 && Counts.size() > 1) {
            --n1;
            Counts.erase(Counts.begin());
        }

        double temp = MultisetCombRowNumFast(n1, r1, Counts);

        for (; temp <= dblIdx; ++j) {
            dblIdx -= temp;
            TempReps[j] = 0;

            if (static_cast<int>(Counts.size()) == (n - j)) {
                --n1;
                Counts.erase(Counts.begin());
            }

            --Counts.front();
            if (Counts.front() == 0 && Counts.size() > 1) {
                --n1;
                Counts.erase(Counts.begin());
            }

            temp = MultisetCombRowNumFast(n1, r1, Counts);
        }

        res[k] = j;
        --TempReps[j];
        if (TempReps[j] <= 0) ++j;
    }

    return res;
}

std::vector<int> nthCombGmp(int n, int m, double dblIdx,
                            const mpz_class &mpzIdx,
                            const std::vector<int> &Reps) {

    mpz_class idx(mpzIdx);
    mpz_class temp;

    std::vector<int> res(m);
    nChooseKGmp(temp, n - 1, m - 1);

    for (int k = 0, j = 0, n1 = n - 1, r1 = m - 1;
         k < m; ++k, --n1, --r1, ++j) {

        for (int s = n1 - r1; cmp(temp, idx) <= 0; --s, ++j, --n1) {
            idx -= temp;
            temp *= s;
            mpz_divexact_ui(temp.get_mpz_t(), temp.get_mpz_t(), n1);
        }

        temp *= r1;
        if (n1 > 0) mpz_divexact_ui(temp.get_mpz_t(), temp.get_mpz_t(), n1);
        res[k] = j;
    }

    return res;
}

std::vector<int> nthCombRepGmp(int n, int m, double dblIdx,
                               const mpz_class &mpzIdx,
                               const std::vector<int> &Reps) {

    mpz_class idx(mpzIdx);
    mpz_class temp;

    std::vector<int> res(m);
    NumCombsWithRepGmp(temp, n, m - 1);

    for (int k = 0, j = 0, n1 = n, r1 = m - 1; k < m; ++k, --r1) {
        for (; cmp(temp, idx) <= 0; ++j, --n1) {
            idx -= temp;
            temp *= (n1 - 1);
            mpz_divexact_ui(temp.get_mpz_t(), temp.get_mpz_t(), n1 + r1 - 1);
        }

        temp *= r1;

        if ((n1 + r1) > 2) {
            mpz_divexact_ui(temp.get_mpz_t(), temp.get_mpz_t(), n1 + r1 - 1);
        }

        res[k] = j;
    }

    return res;
}

// N.B. In the non-Gmp version, Counts is a vector. Here, we use a deque.
// We have empirically determined that using a vector with smaller cases is
// computationally faster than using a deque, even though the algorithmic
// complexity is larger. However, with larger cases, deque proves to be
// more efficient.
std::vector<int> nthCombMultGmp(int n, int m, double dblIdx,
                                const mpz_class &mpzIdx,
                                const std::vector<int> &Reps) {

    mpz_class idx(mpzIdx);
    mpz_class temp;

    std::vector<int> res(m);
    std::deque<int> Counts(Reps.cbegin(), Reps.cend());
    std::vector<int> TempReps(Reps.cbegin(), Reps.cend());

    for (int k = 0, n1 = n, j = 0, r1 = m - 1; k < m; ++k, --r1) {

        --Counts.front();
        if (Counts.size() > 1 && Counts.front() == 0) {
            --n1;
            Counts.pop_front();
        }

        MultisetCombRowNumGmp(temp, n1, r1, Counts);

        for (; cmp(temp, idx) <= 0; ++j) {
            idx -= temp;
            TempReps[j] = 0;

            if (static_cast<int>(Counts.size()) == (n - j)) {
                --n1;
                Counts.pop_front();
            }

            --Counts.front();
            if (Counts.size() > 1 && Counts.front() == 0) {
                --n1;
                Counts.pop_front();
            }

            MultisetCombRowNumGmp(temp, n1, r1, Counts);
        }

        res[k] = j;
        --TempReps[j];
        if (TempReps[j] <= 0) ++j;
    }

    return res;
}

nthCombPtr GetNthCombFunc(bool IsMult, bool IsRep, bool IsGmp) {

    if (IsGmp) {
        if (IsMult) {
            return(nthCombPtr(nthCombMultGmp));
        } else if (IsRep) {
            return(nthCombPtr(nthCombRepGmp));
        } else {
            return(nthCombPtr(nthCombGmp));
        }
    } else {
        if (IsMult) {
            return(nthCombPtr(nthCombMult));
        } else if (IsRep) {
            return(nthCombPtr(nthCombRep));
        } else {
            return(nthCombPtr(nthComb));
        }
    }
}
