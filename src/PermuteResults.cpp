#include "Permutations/PermutationResults.h"

template <typename T>
void PermuteResDistinct(T* mat, const std::vector<T> &v,
                        std::vector<int> &z, int n, int m,
                        int nRows, const funcPtr<T> myFun) {

    std::vector<T> vPass(m);
    auto arrPerm = FromCpp14::make_unique<int[]>(n);

    for (int i = 0; i < n; ++i) {
        arrPerm[i] = z[i];
    }

    if (m == n) {
        // Since we are getting all permutations of v, we know that
        // the result of myFun on v will remain the same for the 5
        // functions defined in ConstraintsUtils.h (i.e. order does
        // not matter for min, max, prod, mean, & sum).
        for (int j = 0; j < m; ++j) {
            vPass[j] = v[arrPerm[j]];
            mat[nRows * j] = vPass[j];
        }

        const auto myRes = myFun(vPass, m);
        mat[nRows * m] = myRes;
        nextFullPerm(arrPerm.get(), n - 1);

        for (int count = 1, numR1 = nRows - 1,
             maxInd = n - 1; count < numR1; ++count) {

            for (int j = 0; j < m; ++j) {
                mat[count +  nRows * j] = v[arrPerm[j]];
            }

            mat[count +  nRows * m] = myRes;
            nextFullPerm(arrPerm.get(), maxInd);
        }
    } else {
        for (int count = 0, numR1 = nRows - 1, lastCol = m - 1,
             maxInd = n - 1; count < numR1; ++count) {

            for (int j = 0; j < m; ++j) {
                vPass[j] = v[arrPerm[j]];
                mat[count +  nRows * j] = vPass[j];
            }

            mat[count +  nRows * m] = myFun(vPass, m);
            nextPartialPerm(arrPerm.get(), lastCol, maxInd);
        }
    }

    // Get last permutation
    for (int j = 0; j < m; ++j) {
        vPass[j] = v[arrPerm[j]];
        mat[nRows - 1L +  nRows * j] = vPass[j];
    }

    mat[nRows - 1L +  nRows * m] = myFun(vPass, m);
}

template <typename T>
void PermuteResDistinct(RcppParallel::RMatrix<T> &mat,
                        const std::vector<T> &v, std::vector<int> &z,
                        int n, int m, int strt, int nRows,
                        const funcPtr<T> myFun) {

    std::vector<T> vPass(m);
    auto arrPerm = FromCpp14::make_unique<int[]>(n);

    for (int i = 0; i < n; ++i) {
        arrPerm[i] = z[i];
    }

    if (m == n) {
        // Since we are getting all permutations of v, we know that
        // the result of myFun on v will remain the same for the 5
        // functions defined in ConstraintsUtils.h (i.e. order does
        // not matter for min, max, prod, mean, & sum).
        for (int j = 0; j < m; ++j) {
            vPass[j] = v[arrPerm[j]];
            mat(strt, j) = vPass[j];
        }

        const auto myRes = myFun(vPass, m);
        mat(strt, m) = myRes;
        nextFullPerm(arrPerm.get(), n - 1);

        for (int count = strt + 1, numR1 = nRows - 1,
             maxInd = n - 1; count < numR1; ++count) {

            for (int j = 0; j < m; ++j) {
                mat(count, j) = v[arrPerm[j]];
            }

            mat(count, m) = myRes;
            nextFullPerm(arrPerm.get(), maxInd);
        }
    } else {
        for (int count = strt, numR1 = nRows - 1, lastCol = m - 1,
             maxInd = n - 1; count < numR1; ++count) {

            for (int j = 0; j < m; ++j) {
                vPass[j] = v[arrPerm[j]];
                mat(count, j) = vPass[j];
            }

            mat(count, m) = myFun(vPass, m);
            nextPartialPerm(arrPerm.get(), lastCol, maxInd);
        }
    }

    // Get last permutation
    for (int j = 0; j < m; ++j) {
        vPass[j] = v[arrPerm[j]];
        mat(nRows - 1L, j) = vPass[j];
    }

    mat(nRows - 1L, m) = myFun(vPass, m);
}

template <typename T>
void PermuteResRep(T* mat, const std::vector<T> &v,
                   std::vector<int> &z, int n, int m,
                   int nRows, const funcPtr<T> myFun) {

    std::vector<T> vPass(m);

    for (int count = 0, maxInd = n - 1,
         lastCol = m - 1; count < nRows; ++count) {

        for (int j = 0; j < m; ++j) {
            vPass[j] = v[z[j]];
            mat[count + nRows * j] = vPass[j];
        }

        mat[count + nRows * m] = myFun(vPass, m);

        for (int i = lastCol; i >= 0; --i) {
            if (z[i] != maxInd) {
                ++z[i];
                break;
            } else {
                z[i] = 0;
            }
        }
    }
}

template <typename T>
void PermuteResRep(RcppParallel::RMatrix<T> &mat,
                   const std::vector<T> &v,
                   std::vector<int> &z, int n, int m,
                   int strt, int nRows, const funcPtr<T> myFun) {

    std::vector<T> vPass(m);

    for (int count = strt, maxInd = n - 1,
         lastCol = m - 1; count < nRows; ++count) {

        for (int j = 0; j < m; ++j) {
            vPass[j] = v[z[j]];
            mat(count, j) = vPass[j];
        }

        mat(count, m) = myFun(vPass, m);

        for (int i = lastCol; i >= 0; --i) {
            if (z[i] != maxInd) {
                ++z[i];
                break;
            } else {
                z[i] = 0;
            }
        }
    }
}

template <typename T>
void MultisetPermRes(T* mat, const std::vector<T> &v,
                     std::vector<int> &z, int n, int m, int nRows,
                     const std::vector<int> &freqs, const funcPtr<T> myFun) {

    const int lenFreqs = freqs.size();
    auto arrPerm = FromCpp14::make_unique<int[]>(lenFreqs);
    std::vector<T> vPass(m);

    for (int j = 0; j < lenFreqs; ++j) {
        arrPerm[j] = z[j];
    }

    if (m == lenFreqs) {
        // Since we are getting all permutations of v, we know that
        // the result of myFun on v will remain the same for the 5
        // functions defined in ConstraintsUtils.h (i.e. order does
        // not matter for min, max, prod, mean, & sum).
        for (int j = 0; j < m; ++j) {
            vPass[j] = v[arrPerm[j]];
            mat[j * nRows] = vPass[j];
        }

        const auto myRes = myFun(vPass, m);
        mat[m * nRows] = myRes;
        nextFullPerm(arrPerm.get(), lenFreqs - 1);

        for (int count = 1, numR1 = nRows - 1,
             maxInd = lenFreqs - 1; count < numR1; ++count) {

            for (int j = 0; j < m; ++j) {
                mat[count + j * nRows] = v[arrPerm[j]];
            }

            mat[count + m * nRows] = myRes;
            nextFullPerm(arrPerm.get(), maxInd);
        }
    } else {
        for (int count = 0, numR1 = nRows - 1, lastCol = m - 1,
             maxInd = lenFreqs - 1; count < numR1; ++count) {

            for (int j = 0; j < m; ++j) {
                vPass[j] = v[arrPerm[j]];
                mat[count + j * nRows] = vPass[j];
            }

            mat[count + m * nRows] = myFun(vPass, m);
            nextPartialPerm(arrPerm.get(), lastCol, maxInd);
        }
    }

    // Get last permutation
    for (int j = 0; j < m; ++j) {
        vPass[j] = v[arrPerm[j]];
        mat[nRows - 1 + j * nRows] = vPass[j];
    }

    mat[nRows - 1 + m * nRows] = myFun(vPass, m);
}

template <typename T>
void MultisetPermRes(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                     std::vector<int> &z, int n, int m, int strt, int nRows,
                     const std::vector<int> &freqs, const funcPtr<T> myFun) {

    const int lenFreqs = freqs.size();
    auto arrPerm = FromCpp14::make_unique<int[]>(lenFreqs);
    std::vector<T> vPass(m);

    for (int j = 0; j < lenFreqs; ++j) {
        arrPerm[j] = z[j];
    }

    if (m == lenFreqs) {
        // Since we are getting all permutations of v, we know that
        // the result of myFun on v will remain the same for the 5
        // functions defined in ConstraintsUtils.h (i.e. order does
        // not matter for min, max, prod, mean, & sum).
        for (int j = 0; j < m; ++j) {
            vPass[j] = v[arrPerm[j]];
            mat(strt, j) = vPass[j];
        }

        const auto myRes = myFun(vPass, m);
        mat(strt, m) = myRes;
        nextFullPerm(arrPerm.get(), lenFreqs - 1);

        for (int count = strt + 1, numR1 = nRows - 1,
             maxInd = lenFreqs - 1; count < numR1; ++count) {

            for (int j = 0; j < m; ++j) {
                mat(count, j) = v[arrPerm[j]];
            }

            mat(count, m) = myRes;
            nextFullPerm(arrPerm.get(), maxInd);
        }
    } else {
        for (int count = strt, numR1 = nRows - 1, lastCol = m - 1,
             maxInd = lenFreqs - 1; count < numR1; ++count) {

            for (int j = 0; j < m; ++j) {
                vPass[j] = v[arrPerm[j]];
                mat(count, j) = vPass[j];
            }

            mat(count, m) = myFun(vPass, m);
            nextPartialPerm(arrPerm.get(), lastCol, maxInd);
        }
    }

    // Get last permutation
    for (int j = 0; j < m; ++j) {
        vPass[j] = v[arrPerm[j]];
        mat(nRows - 1, j) = vPass[j];
    }

    mat(nRows - 1, m) = myFun(vPass, m);
}

template void PermuteResDistinct(int*, const std::vector<int>&,
                                 std::vector<int>&, int, int, int,
                                 const funcPtr<int>);
template void PermuteResDistinct(double*, const std::vector<double>&,
                                 std::vector<int>&, int, int, int,
                                 const funcPtr<double>);
template void PermuteResDistinct(RcppParallel::RMatrix<int>&,
                                 const std::vector<int>&, std::vector<int>&,
                                 int, int, int, int, const funcPtr<int>);
template void PermuteResDistinct(RcppParallel::RMatrix<double>&,
                                 const std::vector<double>&,
                                 std::vector<int>&, int, int, int,
                                 int, const funcPtr<double>);

template void PermuteResRep(int*, const std::vector<int>&, std::vector<int>&,
                            int, int, int, const funcPtr<int>);
template void PermuteResRep(double*, const std::vector<double>&,
                            std::vector<int>&,
                            int, int, int, const funcPtr<double>);
template void PermuteResRep(RcppParallel::RMatrix<int>&,
                            const std::vector<int>&,
                            std::vector<int>&, int, int, int,
                            int, const funcPtr<int>);
template void PermuteResRep(RcppParallel::RMatrix<double>&,
                            const std::vector<double>&,
                            std::vector<int>&, int, int, int,
                            int, const funcPtr<double>);

template void MultisetPermRes(int*, const std::vector<int>&,
                              std::vector<int>&, int, int, int,
                              const std::vector<int>&, const funcPtr<int>);
template void MultisetPermRes(double*, const std::vector<double>&,
                              std::vector<int>&, int, int, int,
                              const std::vector<int>&,
                              const funcPtr<double>);
template void MultisetPermRes(RcppParallel::RMatrix<int>&,
                              const std::vector<int>&,
                              std::vector<int>&, int, int, int, int,
                              const std::vector<int>&, const funcPtr<int>);
template void MultisetPermRes(RcppParallel::RMatrix<double>&,
                              const std::vector<double>&,
                              std::vector<int>&, int, int, int, int,
                              const std::vector<int>&,
                              const funcPtr<double>);
