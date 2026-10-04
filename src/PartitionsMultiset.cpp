#include "Partitions/NextPartition.h"
#include "PopulateUtils.h"
#include "RMatrix.h"

template <typename T>
int PartsGenMultiset(T* mat, const std::vector<T> &v,
                     const std::vector<int> &Reps, std::vector<int> &z,
                     std::size_t width, int lastElem,
                     int lastCol, std::size_t nRows) {

    int b = 0;
    int p = 0;
    int e = 0;

    std::vector<int> rpsCnt(Reps.cbegin(), Reps.cend());
    PrepareMultisetPart(rpsCnt, z, b, p, e, lastCol, lastElem);

    const int lastRow = nRows - 1;

    for (int count = 0; count < lastRow; ++count,
        NextMultisetGenPart(rpsCnt, z, e, b, p, lastCol, lastElem)) {

        for (std::size_t k = 0; k < width; ++k) {
            mat[count + nRows * k] = v[z[k]];
        }
    }

    for (std::size_t k = 0; k < width; ++k) {
        mat[lastRow + nRows * k] = v[z[k]];
    }

    return 1;
}

template <typename T>
int PartsGenMultiset(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                     const std::vector<int> &Reps, std::vector<int> &z,
                     int strt, std::size_t width, int lastElem,
                     int lastCol, std::size_t nRows) {

    int b = 0;
    int p = 0;
    int e = 0;

    std::vector<int> rpsCnt(Reps.cbegin(), Reps.cend());
    PrepareMultisetPart(rpsCnt, z, b, p, e, lastCol, lastElem);

    const int lastRow = nRows - 1;

    for (int count = strt; count < lastRow; ++count,
        NextMultisetGenPart(rpsCnt, z, e, b, p, lastCol, lastElem)) {

        for (std::size_t k = 0; k < width; ++k) {
            mat(count, k) = v[z[k]];
        }
    }

    for (std::size_t k = 0; k < width; ++k) {
        mat[lastRow + nRows * k] = v[z[k]];
    }

    return 1;
}

template <typename T>
int PartsGenPermMultiset(T* mat, const std::vector<T> &v,
                         const std::vector<int> &Reps, std::vector<int> &z,
                         std::size_t width, int lastElem,
                         int lastCol, std::size_t nRows) {

    int b = 0;
    int p = 0;
    int e = 0;

    std::vector<int> rpsCnt(Reps.cbegin(), Reps.cend());
    PrepareMultisetPart(rpsCnt, z, b, p, e, lastCol, lastElem);

    for (std::size_t count = 0;;
         NextMultisetGenPart(rpsCnt, z, e, b, p, lastCol, lastElem)) {

        PopulateMatrix(mat, v, z, count, width, nRows, false);
        if (count >= nRows) {break;}
    }

    return 1;
}

template int PartsGenMultiset(int*, const std::vector<int>&,
                              const std::vector<int>&, std::vector<int>&,
                              std::size_t, int, int, std::size_t);
template int PartsGenMultiset(double*,
                              const std::vector<double>&,
                              const std::vector<int>&, std::vector<int>&,
                              std::size_t, int, int, std::size_t);

template int PartsGenMultiset(RcppParallel::RMatrix<int>&,
                              const std::vector<int>&,
                              const std::vector<int>&, std::vector<int>&,
                              int, std::size_t, int, int, std::size_t);
template int PartsGenMultiset(RcppParallel::RMatrix<double>&,
                              const std::vector<double>&,
                              const std::vector<int>&, std::vector<int>&,
                              int, std::size_t, int, int, std::size_t);

template int PartsGenPermMultiset(
    int*, const std::vector<int>&, const std::vector<int>&,
    std::vector<int>&, std::size_t, int, int, std::size_t
);
template int PartsGenPermMultiset(
    double*, const std::vector<double>&, const std::vector<int>&,
    std::vector<int>&, std::size_t, int, int, std::size_t
);
