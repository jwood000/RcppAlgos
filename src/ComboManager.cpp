#include "Combinations/CombinationsMultiset.h"
#include "Combinations/CombinationsDistinct.h"
#include "Combinations/CombinationsRep.h"
#include "RMatrix.h"

template <typename T>
void ComboManager(T* mat, const std::vector<T> &v,
                  std::vector<int> &z, int n, int m, int nRows,
                  const std::vector<int> &freqs, bool IsMult, bool IsRep) {

    if (IsMult) {
        MultisetCombination(mat, v, z, n, m, nRows, freqs);
    } else if (IsRep) {
        CombinationsRep(mat, v, z, n, m, nRows);
    } else {
        CombinationsDistinct(mat, v, z, n, m, nRows);
    }
}

template <typename T>
void ComboParallel(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                   std::vector<int> &z, int n, int m, int strt, int nRows,
                   const std::vector<int> &freqs, bool IsMult, bool IsRep) {

    if (IsMult) {
        MultisetCombination(mat, v, z, n, m, strt, nRows, freqs);
    } else if (IsRep) {
        CombinationsRep(mat, v, z, n, m, strt, nRows);
    } else {
        CombinationsDistinct(mat, v, z, n, m, strt, nRows);
    }
}

void ComboCharacter(SEXP mat, SEXP v, std::vector<int> &z, int n,
                    int m, int nRows, const std::vector<int> &freqs,
                    bool IsMult, bool IsRep) {

    if (IsMult) {
        MultisetCombination(mat, v, z, n, m, nRows, freqs);
    } else if (IsRep) {
        CombinationsRep(mat, v, z, n, m, nRows);
    } else {
        CombinationsDistinct(mat, v, z, n, m, nRows);
    }
}

template void ComboManager(int*, const std::vector<int>&,
                           std::vector<int>&, int, int, int,
                           const std::vector<int>&, bool, bool);

template void ComboManager(double*, const std::vector<double>&,
                           std::vector<int>&, int, int, int,
                           const std::vector<int>&, bool, bool);

template void ComboManager(Rbyte*, const std::vector<Rbyte>&,
                           std::vector<int>&, int, int, int,
                           const std::vector<int>&, bool, bool);

template void ComboManager(Rcomplex*, const std::vector<Rcomplex>&,
                           std::vector<int>&, int, int, int,
                           const std::vector<int>&, bool, bool);

template void ComboParallel(RcppParallel::RMatrix<int>&, const std::vector<int>&,
                            std::vector<int>&, int, int, int, int,
                            const std::vector<int>&, bool, bool);

template void ComboParallel(RcppParallel::RMatrix<double>&, const std::vector<double>&,
                            std::vector<int>&, int, int, int, int,
                            const std::vector<int>&, bool, bool);
