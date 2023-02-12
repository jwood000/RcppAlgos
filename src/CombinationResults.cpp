#include "Combinations/CombinationResults.h"

template <typename T>
void ComboResDistinct(T* mat, const std::vector<T> &v,
                      std::vector<int> &z, std::size_t n, std::size_t m,
                      std::size_t nRows, const funcPtr<T> myFun) {

    std::vector<T> vPass(m);

    for (std::size_t count = 0, m1 = m - 1, nMinusM = n - m; count < nRows; ) {

        std::size_t numIter = n - z[m1];

        if (numIter + count > nRows)
            numIter = nRows - count;

        for (std::size_t i = 0; i < numIter; ++i, ++count, ++z[m1]) {
            for (std::size_t j = 0; j < m; ++j) {
                vPass[j] = v[z[j]];
                mat[count + nRows * j] = vPass[j];
            }

            mat[count + nRows * m] = myFun(vPass, m);
        }

        nextCombSec(z, m1, nMinusM);
    }
}

template <typename T>
void ComboResDistinct(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                      std::vector<int> &z, std::size_t n, std::size_t m,
                      std::size_t strt, std::size_t nRows,
                      const funcPtr<T> myFun) {

    std::vector<T> vPass(m);

    for (std::size_t count = strt, m1 = m - 1, nMinusM = n - m; count < nRows; ) {

        std::size_t numIter = n - z[m1];

        if (numIter + count > nRows)
            numIter = nRows - count;

        for (std::size_t i = 0; i < numIter; ++i, ++count, ++z[m1]) {
            for (std::size_t j = 0; j < m; ++j) {
                vPass[j] = v[z[j]];
                mat(count, j) = vPass[j];
            }

            mat(count, m) = myFun(vPass, m);
        }

        nextCombSec(z, m1, nMinusM);
    }
}

template <typename T>
void ComboResRep(T* mat, const std::vector<T> &v,
                 std::vector<int> &z, std::size_t n, std::size_t m,
                 std::size_t nRows, const funcPtr<T> myFun) {

    std::vector<T> vPass(m);

    for (std::size_t count = 0, m1 = m - 1, n1 = n - 1; count < nRows; ) {

        std::size_t numIter = n - z[m1];

        if (numIter + count > nRows)
            numIter = nRows - count;

        for (std::size_t i = 0; i < numIter; ++i, ++count, ++z[m1]) {
            for (std::size_t j = 0; j < m; ++j) {
                vPass[j] = v[z[j]];
                mat[count + nRows * j] = vPass[j];
            }

            mat[count + nRows * m] = myFun(vPass, m);
        }

        nextCombSecRep(z, m1, n1);
    }
}

template <typename T>
void ComboResRep(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                 std::vector<int> &z, std::size_t n, std::size_t m,
                 std::size_t strt, std::size_t nRows, const funcPtr<T> myFun) {

    std::vector<T> vPass(m);

    for (std::size_t count = strt, m1 = m - 1, n1 = n - 1; count < nRows; ) {

        std::size_t numIter = n - z[m1];

        if (numIter + count > nRows)
            numIter = nRows - count;

        for (std::size_t i = 0; i < numIter; ++i, ++count, ++z[m1]) {
            for (std::size_t j = 0; j < m; ++j) {
                vPass[j] = v[z[j]];
                mat(count, j) = vPass[j];
            }

            mat(count, m) = myFun(vPass, m);
        }

        nextCombSecRep(z, m1, n1);
    }
}

template <typename T>
void MultisetComboResult(T* mat, const std::vector<T> &v,
                         std::vector<int> &z, std::size_t n, std::size_t m,
                         std::size_t nRows, const std::vector<int> &freqs,
                         const funcPtr<T> myFun) {

    std::vector<int> zIndex(n);
    std::vector<T> vPass(m);

    for (std::size_t i = 0; i < n; ++i) {
        zIndex[i] = std::find(freqs.cbegin(),
                              freqs.cend(), i) - freqs.cbegin();
    }

    // pentExtreme is the location in freqs that represents
    // the maximal value of the second to the last element

    for (std::size_t count = 0, m1 = m - 1,
         pentExtreme = freqs.size() - m; count < nRows;) {

        std::size_t numIter = n - z[m1];

        if (numIter + count > nRows)
            numIter = nRows - count;

        for (std::size_t i = 0; i < numIter; ++i, ++count, ++z[m1]) {
            for (std::size_t j = 0; j < m; ++j) {
                vPass[j] = v[z[j]];
                mat[count + nRows * j] = vPass[j];
            }

            mat[count + nRows * m] = myFun(vPass, m);
        }

        nextCombSecMulti(freqs, zIndex, z, m1, pentExtreme);
    }
}

template <typename T>
void MultisetComboResult(RcppParallel::RMatrix<T> &mat,
                         const std::vector<T> &v, std::vector<int> &z,
                         std::size_t n, std::size_t m, std::size_t strt,
                         std::size_t nRows, const std::vector<int> &freqs,
                         const funcPtr<T> myFun) {

    std::vector<int> zIndex(n);
    std::vector<T> vPass(m);

    for (std::size_t i = 0; i < n; ++i) {
        zIndex[i] = std::find(freqs.cbegin(),
                              freqs.cend(), i) - freqs.cbegin();
    }

    // pentExtreme is the location in freqs that represents
    // the maximal value of the second to the last element

    for (std::size_t count = strt, m1 = m - 1,
         pentExtreme = freqs.size() - m; count < nRows;) {

        std::size_t numIter = n - z[m1];

        if (numIter + count > nRows)
            numIter = nRows - count;

        for (std::size_t i = 0; i < numIter; ++i, ++count, ++z[m1]) {
            for (std::size_t j = 0; j < m; ++j) {
                vPass[j] = v[z[j]];
                mat(count, j) = vPass[j];
            }

            mat(count, m) = myFun(vPass, m);
        }

        nextCombSecMulti(freqs, zIndex, z, m1, pentExtreme);
    }
}

template void ComboResDistinct(int*, const std::vector<int>&,
                               std::vector<int>&, std::size_t, std::size_t,
                               std::size_t, const funcPtr<int>);
template void ComboResDistinct(double*, const std::vector<double>&,
                               std::vector<int>&, std::size_t, std::size_t,
                               std::size_t, const funcPtr<double>);
template void ComboResDistinct(RcppParallel::RMatrix<int>&,
                               const std::vector<int>&, std::vector<int>&,
                               std::size_t, std::size_t, std::size_t,
                               std::size_t, const funcPtr<int>);
template void ComboResDistinct(RcppParallel::RMatrix<double>&,
                               const std::vector<double>&, std::vector<int>&,
                               std::size_t, std::size_t, std::size_t,
                               std::size_t, const funcPtr<double>);

template void ComboResRep(int*, const std::vector<int>&, std::vector<int>&,
                          std::size_t, std::size_t, std::size_t,
                          const funcPtr<int>);
template void ComboResRep(double*, const std::vector<double>&,
                          std::vector<int>&, std::size_t, std::size_t,
                          std::size_t, const funcPtr<double>);
template void ComboResRep(RcppParallel::RMatrix<int>&, const std::vector<int>&,
                          std::vector<int>&, std::size_t, std::size_t,
                          std::size_t, std::size_t, const funcPtr<int>);
template void ComboResRep(RcppParallel::RMatrix<double>&,
                          const std::vector<double>&, std::vector<int>&,
                          std::size_t, std::size_t, std::size_t,
                          std::size_t, const funcPtr<double>);

template void MultisetComboResult(int*, const std::vector<int>&,
                                  std::vector<int>&, std::size_t, std::size_t,
                                  std::size_t, const std::vector<int>&,
                                  const funcPtr<int>);
template void MultisetComboResult(double*, const std::vector<double>&,
                                  std::vector<int>&, std::size_t, std::size_t,
                                  std::size_t, const std::vector<int>&,
                                  const funcPtr<double>);
template void MultisetComboResult(RcppParallel::RMatrix<int>&,
                                  const std::vector<int>&, std::vector<int>&,
                                  std::size_t, std::size_t, std::size_t,
                                  std::size_t, const std::vector<int>&,
                                  const funcPtr<int>);
template void MultisetComboResult(RcppParallel::RMatrix<double>&,
                                  const std::vector<double>&, std::vector<int>&,
                                  std::size_t, std::size_t, std::size_t,
                                  std::size_t, const std::vector<int>&,
                                  const funcPtr<double>);
