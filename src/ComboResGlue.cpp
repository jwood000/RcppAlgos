#include "Combinations/CombinationResults.h"

template <typename T>
void ComboResStd(T* mat, const std::vector<T> &v,
                 std::vector<int> &z, int n, int m,
                 int nRows, bool IsMult, bool IsRep,
                 const std::vector<int> &freqs,
                 const funcPtr<T> myFun) {

    if (IsMult) {
        MultisetComboResult(mat, v, z, n, m, nRows, freqs, myFun);
    } else if (IsRep) {
        ComboResRep(mat, v, z, n, m, nRows, myFun);
    } else {
        ComboResDistinct(mat, v, z, n, m, nRows, myFun);
    }
}

template <typename T>
void ComboResPar(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                 std::vector<int> &z, int n, int m, int strt, int nRows,
                 const std::vector<int> &freqs, const funcPtr<T> myFun,
                 bool IsMult, bool IsRep) {

    if (IsMult) {
        MultisetComboResult(mat, v, z, n, m, strt, nRows, freqs, myFun);
    } else if (IsRep) {
        ComboResRep(mat, v, z, n, m, strt, nRows, myFun);
    } else {
        ComboResDistinct(mat, v, z, n, m, strt, nRows, myFun);
    }
}

template void ComboResPar(RcppParallel::RMatrix<int>&, const std::vector<int>&,
                            std::vector<int>&, int, int, int, int,
                            const std::vector<int>&, const funcPtr<int>, bool, bool);

template void ComboResPar(RcppParallel::RMatrix<double>&, const std::vector<double>&,
                            std::vector<int>&, int, int, int, int,
                            const std::vector<int>&, const funcPtr<double>, bool, bool);

template void ComboResStd(int*, const std::vector<int>&, std::vector<int>&,
                          int, int, int, bool, bool,
                          const std::vector<int>&, const funcPtr<int>);

template void ComboResStd(double*, const std::vector<double>&, std::vector<int>&,
                          int, int, int, bool, bool,
                          const std::vector<int>&, const funcPtr<double>);
