#include "Permutations/PermuteHelper.h"
#include <cmath>

constexpr std::size_t unrollSize = 8;

template <typename T>
void RepUnroller(T* mat, T val, std::size_t strt,
                 std::size_t last, std::size_t lastUnroll) {

    for (std::size_t i = strt; i < lastUnroll; i += unrollSize) {
        mat[i] = val;
        mat[i + 1] = val;
        mat[i + 2] = val;
        mat[i + 3] = val;
        mat[i + 4] = val;
        mat[i + 5] = val;
        mat[i + 6] = val;
        mat[i + 7] = val;
    }

    for (std::size_t i = lastUnroll; i < last; ++i) {
        mat[i] = val;
    }
}

template <typename T>
void StandardUnroller(T* mat, const int *const indexMat,
                      const std::vector<T> &v, std::size_t m, std::size_t strt,
                      std::size_t last, std::size_t first,
                      std::size_t lastUnroll, std::size_t nRows) {

    for (std::size_t j = first * nRows, k = 0; j < (m * nRows); j += nRows) {
        for (std::size_t i = strt; i < lastUnroll;
             i += unrollSize, k += unrollSize) {
            mat[i + j] = v[indexMat[k]];
            mat[i + 1 + j] = v[indexMat[k + 1]];
            mat[i + 2 + j] = v[indexMat[k + 2]];
            mat[i + 3 + j] = v[indexMat[k + 3]];
            mat[i + 4 + j] = v[indexMat[k + 4]];
            mat[i + 5 + j] = v[indexMat[k + 5]];
            mat[i + 6 + j] = v[indexMat[k + 6]];
            mat[i + 7 + j] = v[indexMat[k + 7]];
        }

        for (std::size_t i = lastUnroll; i < last; ++i, ++k) {
            mat[i + j] = v[indexMat[k]];
        }
    }
}

template <typename T>
void PermuteWorker(T* mat, const int *const indexMat,
                   const std::vector<T> &v, int m,
                   int strt, int last, int ind, int first,
                   int unrollRem, bool IsRep, int nRows) {

    const int lastUnroll = last - unrollRem;

    // For IsRep case, we are not setting the first column because we
    // know that it simply increments. This is taken into account in
    // the indexMat preparation.
    if (IsRep) {
        RepUnroller(mat, v[ind], strt, last, lastUnroll);
    }

    StandardUnroller(mat, indexMat, v, m, strt,
                     last, first, lastUnroll, nRows);
}

template <typename T>
void PermuteLoadIndex(T* mat, int *const indexMat,
                      const std::vector<T> &v, std::vector<int> &z,
                      std::size_t n, std::size_t m, std::size_t segment,
                      bool IsRep, std::size_t nRows) {

    if (IsRep) {
        for (std::size_t count = 0, maxInd = n - 1,
             lastCol = m - 1; count < segment; ++count) {

            // N.B. In PermuteGeneral we start j at 0
            for (std::size_t j = 1, k = count; j < m; ++j, k += segment) {
                mat[count + j * nRows] = v[z[j]];
                indexMat[k] = z[j];
            }

            mat[count] = v[z.front()];

            // N.B. In PermuteGeneral we decrement i until i == 0
            for (int i = lastCol; i > 0; --i) {
                if (z[i] != static_cast<int>(maxInd)) {
                    ++z[i];
                    break;
                } else {
                    z[i] = 0;
                }
            }
        }
    } else {
        auto arrPerm = std::make_unique<int[]>(n);

        for (std::size_t i = 0; i < n; ++i) {
            arrPerm[i] = z[i];
        }

        if (m == n) {
            for (std::size_t count = 0, maxInd = n - 1;
                 count < segment; ++count) {
                for (std::size_t j = 0, k = count; j < m; ++j, k += segment) {
                    mat[count + j * nRows] = v[arrPerm[j]];
                    indexMat[k] = arrPerm[j];
                }

                nextFullPerm(arrPerm.get(), maxInd);
            }
        } else {
            for (std::size_t count = 0, maxInd = n - 1,
                 lastCol = m - 1; count < segment; ++count) {

                for (std::size_t j = 0, k = count; j < m; ++j, k += segment) {
                    mat[count + j * nRows] = v[arrPerm[j]];
                    indexMat[k] = arrPerm[j];
                }

                nextPartialPerm(arrPerm.get(), lastCol, maxInd);
            }
        }
    }
}

template <typename T>
void PermuteOptimized(T* mat, const std::vector<T> &v, std::vector<int> &z,
                      std::size_t n, std::size_t m, std::size_t nRows,
                      bool IsRep) {

    const std::size_t first = (IsRep) ? 1 : 0;
    const std::size_t segment = (IsRep) ? std::pow(static_cast<double>(n),
                                                   static_cast<double>(m - 1)) :
                        NumPermsNoRep(n - 1, m - 1);

    const std::size_t indexMatSize = segment * (m - first);
    auto indexMat = std::make_unique<int[]>(indexMatSize);
    PermuteLoadIndex(mat, indexMat.get(), v, z, n, m, segment, IsRep, nRows);

    std::size_t ind = 1;
    std::size_t strt = segment;
    std::size_t last = strt + segment;
    std::size_t unrollRem = segment % unrollSize;
    std::vector<T> vCopy(v.cbegin(), v.cend());

    for (; last <= nRows; strt += segment, last += segment, ++ind) {
        if (!IsRep) {
            std::swap(vCopy.front(), vCopy[ind]);
        }

        PermuteWorker(mat, indexMat.get(), vCopy, m, strt,
                      last, ind, first, unrollRem, IsRep, nRows);
    }

    if (ind < vCopy.size() && strt < nRows) {
        const std::size_t skip = last - nRows;
        unrollRem = nRows % unrollSize;

        if (IsRep) {
            const std::size_t lastUnroll = nRows - unrollRem;
            RepUnroller(mat, vCopy[ind], strt, nRows, lastUnroll);
        } else {
            std::swap(vCopy.front(), vCopy[ind]);
        }

        for (std::size_t j = first * nRows, k = 0; j < (m * nRows);
             j += nRows, k += skip) {

            for (std::size_t i = strt; i < nRows; ++i, ++k) {
                mat[i + j] = vCopy[indexMat[k]];
            }
        }
    }
}

template void PermuteOptimized(int*, const std::vector<int>&,
                               std::vector<int>&, std::size_t, std::size_t,
                               std::size_t, bool);
template void PermuteOptimized(double*, const std::vector<double>&,
                               std::vector<int>&, std::size_t, std::size_t,
                               std::size_t, bool);
template void PermuteOptimized(Rbyte*, const std::vector<Rbyte>&,
                               std::vector<int>&, std::size_t, std::size_t,
                               std::size_t, bool);
template void PermuteOptimized(Rcomplex*, const std::vector<Rcomplex>&,
                               std::vector<int>&, std::size_t, std::size_t,
                               std::size_t, bool);
