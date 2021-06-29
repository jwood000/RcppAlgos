#include "Permutations/PermutationsMultiset.h"
#include "Permutations/PermutationsDistinct.h"
#include "Permutations/PermuteOptimized.h"
#include "Permutations/PermutationsRep.h"

template <typename T>
void PermuteManager(T* mat, const std::vector<T> &v,
                    std::vector<int> &z, int n, int m, int nRows,
                    int phaseOne, bool generalRet, bool IsMult,
                    bool IsRep, const std::vector<int> &freqs) {

    if (generalRet) {
        if (IsMult) {
            PermuteMultiset(mat, v, z, n, m, nRows, freqs);
        } else if (IsRep) {
            PermuteRep(mat, v, z, n, m, nRows);
        } else {
            PermuteDistinct(mat, v, z, n, m, nRows);
        }
    } else {
        PermuteOptimized(mat, v, z, n, m, nRows, IsRep);
    }
}

template <typename T>
void PermuteParallel(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                     std::vector<int> &z, int n, int m, int strt, int nRows,
                     const std::vector<int> &freqs, bool IsMult, bool IsRep) {

    if (IsMult) {
        PermuteMultiset(mat, v, z, n, m, strt, nRows, freqs);
    } else if (IsRep) {
        PermuteRep(mat, v, z, n, m, strt, nRows);
    } else {
        PermuteDistinct(mat, v, z, n, m, strt, nRows);
    }
}

void PermuteCharacter(SEXP mat, SEXP v, std::vector<int> &z, int n,
                      int m, int nRows, const std::vector<int> &freqs,
                      bool IsMult, bool IsRep) {

    if (IsMult) {
        PermuteMultiset(mat, v, z, n, m, nRows, freqs);
    } else if (IsRep) {
        PermuteRep(mat, v, z, n, m, nRows);
    } else {
        PermuteDistinct(mat, v, z, n, m, nRows);
    }
}

template void PermuteParallel(RcppParallel::RMatrix<int>&,
                              const std::vector<int>&,
                              std::vector<int>&, int, int, int, int,
                              const std::vector<int>&, bool, bool);

template void PermuteParallel(RcppParallel::RMatrix<double>&,
                              const std::vector<double>&,
                              std::vector<int>&, int, int, int, int,
                              const std::vector<int>&, bool, bool);

template void PermuteManager(int*, const std::vector<int>&,
                             std::vector<int>&, int, int, int, int,
                             bool, bool, bool, const std::vector<int>&);

template void PermuteManager(double*, const std::vector<double>&,
                             std::vector<int>&, int, int, int, int,
                             bool, bool, bool, const std::vector<int>&);

template void PermuteManager(Rbyte*, const std::vector<Rbyte>&,
                             std::vector<int>&, int, int, int, int,
                             bool, bool, bool, const std::vector<int>&);

template void PermuteManager(Rcomplex*, const std::vector<Rcomplex>&,
                             std::vector<int>&, int, int, int, int,
                             bool, bool, bool, const std::vector<int>&);
