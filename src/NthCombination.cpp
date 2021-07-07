#include "Combinations/BigComboCount.h"
#include "Combinations/ComboCount.h"

using nthCombPtr = std::vector<int> (*const)(int n, int r, double dblIdx,
                                     mpz_t mpzIdx, const std::vector<int> &Reps);

std::vector<int> nthComb(int n, int r, double dblIdx,
                         mpz_t mpzIdx, const std::vector<int> &Reps) {

    double index1 = dblIdx;
    double index2 = dblIdx;

    std::vector<int> res(r);
    double temp = nChooseK(n - 1, r - 1);

    for (int k = 0, j = 0, n1 = n - 1,
         r1 = r - 1; k < r; ++k, --n1, --r1, ++j) {
        double test = temp;

        for (int s = n1 - r1; test <= index1;
                    --n1, ++j, --s, test += temp) {
            index2 -= temp;
            temp *= s;
            temp /= n1;
        }

        temp *= r1;
        temp /= n1;
        res[k] = j;
        index1 = index2;
    }

    return res;
}

std::vector<int> nthCombRep(int n, int r, double dblIdx,
                            mpz_t mpzIdx, const std::vector<int> &Reps) {

    double index1 = dblIdx;
    double index2 = dblIdx;
    std::vector<int> res(r);
    double temp = NumCombsWithRep(n, r - 1);

    for (int k = 0, j = 0, n1 = n, r1 = r - 1; k < r; ++k, --r1) {
        double test = temp;

        for (; test <= index1; --n1, ++j, test += temp) {
            index2 -= temp;
            temp *= (n1 - 1);
            temp /= (n1 + r1 - 1);
        }

        temp *= r1;
        temp /= (n1 + r1 - 1);
        res[k] = j;
        index1 = index2;
    }

    return res;
}

std::vector<int> nthCombMult(int n, int r, double dblIdx,
                             mpz_t mpzIdx, const std::vector<int> &Reps) {

    double index1 = dblIdx;
    double index2 = dblIdx;
    std::vector<int> res(r);
    std::vector<int> Counts = Reps;
    std::vector<int> TempReps = Reps;

    for (int k = 0, j = 0, n1 = n, r1 = r - 1; k < r; ++k, --r1) {

        --Counts[0];
        if (Counts[0] == 0 && Counts.size() > 1) {
            --n1;
            Counts.erase(Counts.begin());
        }

        double test = MultisetCombRowNumFast(n1, r1, Counts);
        double temp = test;

        for (; test <= index1; ++j, test += temp) {
            index2 -= temp;
            TempReps[j] = 0;

            if (static_cast<int>(Counts.size()) == (n - j)) {
                --n1;
                Counts.erase(Counts.begin());
            }

            --Counts[0];
            if (Counts[0] == 0 && Counts.size() > 1) {
                --n1;
                Counts.erase(Counts.begin());
            }

            temp = MultisetCombRowNumFast(n1, r1, Counts);
        }

        res[k] = j;
        index1 = index2;

        --TempReps[j];
        if (TempReps[j] <= 0) ++j;
    }

    return res;
}

std::vector<int> nthCombGmp(int n, int r, double dblIdx,
                            mpz_t mpzIdx, const std::vector<int> &Reps) {

    mpz_t test;
    mpz_t temp;
    mpz_t index1;
    mpz_t index2;

    mpz_init(test);
    mpz_init(temp);
    mpz_init(index1);
    mpz_init(index2);
    
    mpz_set(index1, mpzIdx);
    mpz_set(index2, mpzIdx);

    std::vector<int> res(r);
    nChooseKGmp(temp, n - 1, r - 1);

    for (int k = 0, j = 0, n1 = n - 1, r1 = r - 1;
         k < r; ++k, --n1, --r1, ++j) {
        mpz_set(test, temp);

        for (int s = n1 - r1; mpz_cmp(test, index1) <= 0; --s, ++j, --n1) {
            mpz_sub(index2, index2, temp);
            mpz_mul_ui(temp, temp, s);
            mpz_divexact_ui(temp, temp, n1);
            mpz_add(test, test, temp);
        }

        mpz_mul_ui(temp, temp, r1);
        if (n1 > 0) mpz_divexact_ui(temp, temp, n1);
        res[k] = j;
        mpz_set(index1, index2);
    }

    mpz_clear(index1); mpz_clear(index2);
    mpz_clear(temp); mpz_clear(test);
    return res;
}

std::vector<int> nthCombRepGmp(int n, int r, double dblIdx,
                               mpz_t mpzIdx, const std::vector<int> &Reps) {

    mpz_t test;
    mpz_t temp;
    mpz_t index1;
    mpz_t index2;

    mpz_init(test);
    mpz_init(temp);
    mpz_init(index1);
    mpz_init(index2);
    
    mpz_set(index1, mpzIdx);
    mpz_set(index2, mpzIdx);

    std::vector<int> res(r);
    NumCombsWithRepGmp(temp, n, r - 1);

    for (int k = 0, j = 0, n1 = n, r1 = r - 1; k < r; ++k, --r1) {
        mpz_set(test, temp);

        for (; mpz_cmp(test, index1) <= 0; ++j, --n1) {
            mpz_sub(index2, index2, temp);
            mpz_mul_ui(temp, temp, n1 - 1);
            mpz_divexact_ui(temp, temp, n1 + r1 - 1);
            mpz_add(test, test, temp);
        }

        mpz_mul_ui(temp, temp, r1);
        if ((n1 + r1) > 1) mpz_divexact_ui(temp, temp, n1 + r1 - 1);
        res[k] = j;
        mpz_set(index1, index2);
    }

    mpz_clear(index1); mpz_clear(index2);
    mpz_clear(temp); mpz_clear(test);
    return res;
}

std::vector<int> nthCombMultGmp(int n, int r, double dblIdx,
                                mpz_t mpzIdx, const std::vector<int> &Reps) {

    mpz_t test;
    mpz_t temp;
    mpz_t index1;
    mpz_t index2;

    mpz_init(test);
    mpz_init(temp);
    mpz_init(index1);
    mpz_init(index2);
    
    mpz_set(index1, mpzIdx);
    mpz_set(index2, mpzIdx);

    std::vector<int> res(r);
    std::vector<int> Counts = Reps;
    std::vector<int> TempReps = Reps;

    for (int k = 0, n1 = n, j = 0, r1 = r - 1; k < r; ++k, --r1) {

        --Counts[0];
        if (Counts[0] == 0 && Counts.size() > 1) {
            --n1;
            Counts.erase(Counts.begin());
        }

        MultisetCombRowNumGmp(temp, n1, r1, Counts);
        mpz_set(test, temp);

        for (; mpz_cmp(test, index1) <= 0; ++j) {
            mpz_sub(index2, index2, temp);
            TempReps[j] = 0;

            if (static_cast<int>(Counts.size()) == (n - j)) {
                --n1;
                Counts.erase(Counts.begin());
            }

            --Counts[0];
            if (Counts[0] == 0 && Counts.size() > 1) {
                --n1;
                Counts.erase(Counts.begin());
            }

            MultisetCombRowNumGmp(temp, n1, r1, Counts);
            mpz_add(test, test, temp);
        }

        res[k] = j;
        mpz_set(index1, index2);

        --TempReps[j];
        if (TempReps[j] <= 0) ++j;
    }

    mpz_clear(index1); mpz_clear(index2);
    mpz_clear(temp); mpz_clear(test);
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
