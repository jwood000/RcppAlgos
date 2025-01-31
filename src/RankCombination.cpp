#include "Combinations/BigComboCount.h"
#include "Combinations/ComboCount.h"

using rankCombPtr = void (*const)(std::vector<int>::iterator iter, int n,
                          int m, double &dblIdx, mpz_class &mpzIdx,
                          const std::vector<int> &Reps);

void rankComb(std::vector<int>::iterator iter, int n, int m,
              double &dblIdx, mpz_class &mpzIdx,
              const std::vector<int> &Reps) {

    dblIdx = 0;
    double temp = nChooseK(n - 1, m - 1);

    for (int k = 0, j = 0, n1 = n - 1, r1 = m - 1;
         k < m; ++k, --n1, --r1, ++j, ++iter) {

        for (int s = n1 - r1, idx = *iter; j < idx; --n1, ++j, --s) {
            dblIdx += temp;
            temp *= s;
            temp /= n1;
        }

        temp *= r1;
        temp /= n1;
    }
}

void rankCombRep(std::vector<int>::iterator iter, int n, int m,
                 double &dblIdx, mpz_class &mpzIdx,
                 const std::vector<int> &Reps) {

    dblIdx = 0;
    double temp = NumCombsWithRep(n, m - 1);

    for (int k = 0, j = 0, n1 = n, r1 = m - 1; k < m; ++k, --r1, ++iter) {
        for (int idx = *iter; j < idx; --n1, ++j) {
            dblIdx += temp;
            temp *= (n1 - 1);
            temp /= (n1 + r1 - 1);
        }

        temp *= r1;
        temp /= (n1 + r1 - 1);
    }
}

void rankCombMult(std::vector<int>::iterator iter, int n, int m,
                  double &dblIdx, mpz_class &mpzIdx,
                  const std::vector<int> &Reps) {

    dblIdx = 0;
    std::vector<int> Counts = Reps;
    std::vector<int> TempReps = Reps;

    for (int k = 0, j = 0, n1 = n, r1 = m - 1; k < m; ++k, --r1, ++iter) {

        ManageCountsVector(Counts, n1);
        double temp = MultisetCombRowNumFast(n1, r1, Counts);

        for (int idx = *iter; j < idx; ++j) {
            dblIdx += temp;
            TempReps[j] = 0;

            if (static_cast<int>(Counts.size()) == (n - j)) {
                --n1;
                Counts.erase(Counts.begin());
            }

            ManageCountsVector(Counts, n1);
            temp = MultisetCombRowNumFast(n1, r1, Counts);
        }

        --TempReps[j];
        if (TempReps[j] <= 0) ++j;
    }
}

void rankCombGmp(std::vector<int>::iterator iter, int n, int m,
                 double &dblIdx, mpz_class &mpzIdx,
                 const std::vector<int> &Reps) {

    mpz_class temp;
    mpzIdx = 0;
    nChooseKGmp(temp, n - 1, m - 1);

    for (int k = 0, j = 0, n1 = n - 1, r1 = m - 1;
         k < m; ++k, --n1, --r1, ++j, ++iter) {

        for (int s = n1 - r1, idx = *iter; j < idx; --s, ++j, --n1) {
            mpzIdx += temp;
            temp *= s;
            mpz_divexact_ui(temp.get_mpz_t(), temp.get_mpz_t(), n1);
        }

        temp *= r1;
        if (n1 > 0) mpz_divexact_ui(temp.get_mpz_t(), temp.get_mpz_t(), n1);
    }
}

void rankCombRepGmp(std::vector<int>::iterator iter, int n, int m,
                    double &dblIdx, mpz_class &mpzIdx,
                    const std::vector<int> &Reps) {

    mpz_class temp;
    mpzIdx = 0;
    NumCombsWithRepGmp(temp, n, m - 1);

    for (int k = 0, j = 0, n1 = n, r1 = m - 1; k < m; ++k, --r1, ++iter) {
        for (int idx = *iter; j < idx; ++j, --n1) {
            mpzIdx += temp;
            temp *= (n1 - 1);
            mpz_divexact_ui(temp.get_mpz_t(), temp.get_mpz_t(), n1 + r1 - 1);
        }

        temp *= r1;
        if ((n1 + r1) > 2) mpz_divexact_ui(temp.get_mpz_t(), temp.get_mpz_t(), n1 + r1 - 1);
    }
}

// N.B. In the non-Gmp version, Counts is a vector. Here, we use a deque.
// We have empirically determined that using a vector with smaller cases is
// computationally faster than using a deque, even though the algorithmic
// complexity is larger. However, with larger cases, deque proves to be
// more efficient.
void rankCombMultGmp(std::vector<int>::iterator iter, int n, int m,
                     double &dblIdx, mpz_class &mpzIdx,
                     const std::vector<int> &Reps) {

    mpz_class temp;
    mpzIdx = 0;

    std::deque<int> Counts(Reps.cbegin(), Reps.cend());
    std::vector<int> TempReps(Reps.cbegin(), Reps.cend());

    for (int k = 0, n1 = n, j = 0, r1 = m - 1; k < m; ++k, --r1, ++iter) {

        ManageCountsDeque(Counts, n1);
        MultisetCombRowNumGmp(temp, n1, r1, Counts);

        for (int idx = *iter; j < idx; ++j) {
            mpzIdx += temp;
            TempReps[j] = 0;

            if (static_cast<int>(Counts.size()) == (n - j)) {
                --n1;
                Counts.pop_front();
            }

            ManageCountsDeque(Counts, n1);
            MultisetCombRowNumGmp(temp, n1, r1, Counts);
        }

        --TempReps[j];
        if (TempReps[j] <= 0) ++j;
    }
}

rankCombPtr GetRankCombFunc(bool IsMult, bool IsRep, bool IsGmp) {

    if (IsGmp) {
        if (IsMult) {
            return(rankCombPtr(rankCombMultGmp));
        } else if (IsRep) {
            return(rankCombPtr(rankCombRepGmp));
        } else {
            return(rankCombPtr(rankCombGmp));
        }
    } else {
        if (IsMult) {
            return(rankCombPtr(rankCombMult));
        } else if (IsRep) {
            return(rankCombPtr(rankCombRep));
        } else {
            return(rankCombPtr(rankComb));
        }
    }
}
