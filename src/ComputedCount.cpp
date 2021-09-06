#include "Permutations/BigPermuteCount.h"
#include "Combinations/BigComboCount.h"
#include "Permutations/PermuteCount.h"
#include "Combinations/ComboCount.h"
#include "CleanConvert.h"
#include <cmath>

double GetComputedRows(bool IsMult, bool IsComb, bool IsRep,
                       int n, int m, const SEXP &Rm,
                       const std::vector<int> &freqs,
                       const std::vector<int> &Reps) {

    double computedRows = 0;

    if (IsMult) {
        if (IsComb) {
            computedRows = MultisetCombRowNum(n, m, Reps);
        } else {
            if (Rf_isNull(Rm) || m == static_cast<int>(freqs.size())) {
                computedRows = NumPermsWithRep(freqs);
            } else {
                computedRows = MultisetPermRowNum(n, m, Reps);
            }
        }
    } else {
        if (IsRep) {
            if (IsComb) {
                computedRows = NumCombsWithRep(n, m);
            } else {
                computedRows = std::pow(static_cast<double>(n),
                                        static_cast<double>(m));
            }
        } else {
            if (IsComb) {
                computedRows = nChooseK(n, m);
            } else {
                computedRows = NumPermsNoRep(n, m);
            }
        }
    }

    return computedRows;
}

void GetComputedRowMpz(mpz_t computedRowsMpz, bool IsMult, bool IsComb,
                       bool IsRep, int n, int m, const SEXP &Rm,
                       const std::vector<int> &freqs,
                       const std::vector<int> &myReps) {

    if (IsMult) {
        if (IsComb) {
            std::deque<int> deqRes(myReps.cbegin(), myReps.cend());
            MultisetCombRowNumGmp(computedRowsMpz, n, m, deqRes);
        } else {
            if (Rf_isNull(Rm) || m == static_cast<int>(freqs.size())) {
                NumPermsWithRepGmp(computedRowsMpz, freqs);
            } else {
                MultisetPermRowNumGmp(computedRowsMpz, n, m, myReps);
            }
        }
    } else {
        if (IsRep) {
            if (IsComb) {
                NumCombsWithRepGmp(computedRowsMpz, n, m);
            } else {
                mpz_ui_pow_ui(computedRowsMpz, n, m);
            }
        } else {
            if (IsComb) {
                nChooseKGmp(computedRowsMpz, n, m);
            } else {
                NumPermsNoRepGmp(computedRowsMpz, n, m);
            }
        }
    }
}
