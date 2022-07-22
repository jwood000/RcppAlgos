#include "Partitions/NextComposition.h"
#include "RMatrix.h"
#include <algorithm>  // std::next_permutation

template <typename T>
void CompsGenRep(T* mat, const std::vector<T> &v,
                 std::vector<int> &z, int width, int nRows) {

    for (int count = 0, lastCol = width - 1; count < nRows; ++count,
         NextCompositionRep<0>(z, lastCol)) {

        for (int k = 0; k < width; ++k) {
            mat[count + nRows * k] = v[z[k]];
        }
    }
}

template <typename T>
void CompsGenRep(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                 std::vector<int> &z, int strt, int width, int nRows) {

    for (int count = strt, lastCol = width - 1; count < nRows; ++count,
         NextCompositionRep<0>(z, lastCol)) {

        for (int k = 0; k < width; ++k) {
            mat(count, k) = v[z[k]];
        }
    }
}

void CompsRep(int* mat, std::vector<int> &z,
              int width, int nRows) {

    for (int count = 0, lastCol = width - 1; count < nRows; ++count,
         NextCompositionRep<1>(z, lastCol)) {

        for (int k = 0; k < width; ++k) {
            mat[count + nRows * k] = z[k];
        }
    }
}

void CompsRep(RcppParallel::RMatrix<int> &mat, std::vector<int> &z,
              int strt, int width, int nRows) {

    for (int count = strt, lastCol = width - 1; count < nRows; ++count,
         NextCompositionRep<1>(z, lastCol)) {

        for (int k = 0; k < width; ++k) {
            mat(count, k) = z[k];
        }
    }
}

template void CompsGenRep(int*, const std::vector<int>&,
                          std::vector<int>&, int, int);
template void CompsGenRep(double*, const std::vector<double>&,
                          std::vector<int>&, int, int);

template void CompsGenRep(RcppParallel::RMatrix<int>&,
                          const std::vector<int>&,
                          std::vector<int>&, int, int, int);
template void CompsGenRep(RcppParallel::RMatrix<double>&,
                          const std::vector<double>&,
                          std::vector<int>&, int, int, int);
