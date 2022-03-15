#include "Permutations/PermutationResults.h"

template <typename T>
void PermuteResStd(T* mat, const std::vector<T> &v,
                   std::vector<int> &z, int n, int m,
                   int nRows, bool IsMult, bool IsRep,
                   const std::vector<int> &freqs,
                   const funcPtr<T> myFun) {

    if (IsMult) {
        MultisetPermRes(mat, v, z, n, m, nRows, freqs, myFun);
    } else if (IsRep) {
        PermuteResRep(mat, v, z, n, m, nRows, myFun);
    } else {
        PermuteResDistinct(mat, v, z, n, m, nRows, myFun);
    }
}

template <typename T>
void PermuteResPar(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                   std::vector<int> &z, int n, int m, int strt, int nRows,
                   const std::vector<int> &freqs, const funcPtr<T> myFun,
                   bool IsMult, bool IsRep) {

    if (IsMult) {
        MultisetPermRes(mat, v, z, n, m, strt, nRows, freqs, myFun);
    } else if (IsRep) {
        PermuteResRep(mat, v, z, n, m, strt, nRows, myFun);
    } else {
        PermuteResDistinct(mat, v, z, n, m, strt, nRows, myFun);
    }
}

template void PermuteResPar(RcppParallel::RMatrix<int>&,
                            const std::vector<int>&,
                            std::vector<int>&, int, int, int, int,
                            const std::vector<int>&, const funcPtr<int>,
                            bool, bool);

template void PermuteResPar(RcppParallel::RMatrix<double>&,
                            const std::vector<double>&,
                            std::vector<int>&, int, int, int, int,
                            const std::vector<int>&, const funcPtr<double>,
                            bool, bool);

template void PermuteResStd(int*, const std::vector<int>&, std::vector<int>&,
                            int, int, int, bool, bool,
                            const std::vector<int>&, const funcPtr<int>);

template void PermuteResStd(double*, const std::vector<double>&,
                            std::vector<int>&, int, int, int, bool, bool,
                            const std::vector<int>&, const funcPtr<double>);
