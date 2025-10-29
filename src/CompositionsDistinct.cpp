#include "Partitions/PartitionsCountDistinct.h"
#include "Partitions/NextComposition.h"
#include "RMatrix.h"
#include <algorithm> // std::set_difference
#include <numeric>

void CompsDistinctWorker(
    int* mat, std::vector<int> &z, std::vector<int> &complement,
    int i1, int i2, int myMax, int tar, std::size_t strt,
    std::size_t width, std::size_t nRows, std::size_t totalRows
) {

    // If the size of complement is 0 then we know we are simply taking
    // permutations of z. Same thing with size 1 but a little reasoning is
    // needed. In NextCompositionDistinct, we are swapping the last two
    // elements in z with elements in complement. If complement only has one
    // element, we are reduced to the case of taking permutations of z.
    if (complement.size() < 2) {
        std::size_t count = strt;

        do {
            for (std::size_t k = 0; k < width; ++k) {
                mat[count + totalRows * k] = z[k];
            }

            ++count;
        } while (std::next_permutation(z.begin(), z.end()) && count < nRows);

        return;
    }

    std::vector<int> idx;
    std::vector<int> tailSum;

    for (std::size_t count = strt, m = width - 1, q = complement.size() - 1,
         lastRow = nRows - 1; count < lastRow; ++count) {

        for (std::size_t k = 0; k < width; ++k) {
            mat[count + totalRows * k] = z[k];
        }

        NextCompositionDistinct(
            z, complement, idx, tailSum, i1, i2, myMax, m, q, tar
        );
    }

    for (std::size_t k = 0; k < width; ++k) {
        mat[nRows - 1 + totalRows * k] = z[k];
    }
}

template <typename T>
void CompsDistinctWorker(
    T* mat, const std::vector<T> &v, std::vector<int> &z,
    std::vector<int> &complement, int i1, int i2, int myMax,
    int tar, std::size_t strt, std::size_t width,
    std::size_t nRows, std::size_t totalRows
) {

    // If the size of complement is 0 then we know we are simply taking
    // permutations of z. Same thing with size 1 but a little reasoning is
    // needed. In NextCompositionDistinct, we are swapping the last two
    // elements in z with elements in complement. If complement only has one
    // element, we are reduced to the case of taking permutations of z.
    if (complement.size() < 2) {
        std::size_t count = strt;

        do {
            for (std::size_t k = 0; k < width; ++k) {
                mat[count + totalRows * k] = v[z[k]];
            }

            ++count;
        } while (std::next_permutation(z.begin(), z.end()) && count < nRows);

        return;
    }

    std::vector<int> idx;
    std::vector<int> tailSum;

    for (std::size_t count = strt, m = width - 1, q = complement.size() - 1,
         lastRow = nRows - 1; count < lastRow; ++count) {

        for (std::size_t k = 0; k < width; ++k) {
            mat[count + totalRows * k] = v[z[k]];
        }

        NextCompositionDistinct(
            z, complement, idx, tailSum, i1, i2, myMax, m, q, tar
        );
    }

    for (std::size_t k = 0; k < width; ++k) {
        mat[nRows - 1 + totalRows * k] = v[z[k]];
    }
}

template <typename T>
void CompsDistinctWorker(
    RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
    std::vector<int> &z, std::vector<int> &complement,
    int i1, int i2, int myMax, int tar, std::size_t strt,
    std::size_t width, std::size_t nRows
) {

    std::vector<int> idx;
    std::vector<int> tailSum;

    for (std::size_t count = strt, m = width - 1, q = complement.size() - 1,
         lastRow = nRows - 1; count < lastRow; ++count) {

        for (std::size_t k = 0; k < width; ++k) {
            mat(count, k) = v[z[k]];
        }

        NextCompositionDistinct(
            z, complement, idx, tailSum, i1, i2, myMax, m, q, tar
        );
    }

    for (std::size_t k = 0; k < width; ++k) {
        mat(nRows - 1, k) = v[z[k]];
    }
}

void CompsDistinctWorker(
    RcppParallel::RMatrix<int> &mat, std::vector<int> &z,
    std::vector<int> &complement, int i1, int i2, int myMax,
    int tar, std::size_t strt, std::size_t width, std::size_t nRows
) {

    std::vector<int> idx;
    std::vector<int> tailSum;

    for (std::size_t count = strt, m = width - 1, q = complement.size() - 1,
         lastRow = nRows - 1; count < lastRow; ++count) {

        for (std::size_t k = 0; k < width; ++k) {
            mat(count, k) = z[k];
        }

        NextCompositionDistinct(
            z, complement, idx, tailSum, i1, i2, myMax, m, q, tar
        );
    }

    for (std::size_t k = 0; k < width; ++k) {
        mat(nRows - 1, k) = z[k];
    }
}

template <typename T>
int CompsGenDistinct(
    T* mat, const std::vector<T> &v, std::vector<int> &z,
    std::size_t width, std::size_t nRows
) {

    int nz;
    int tar;
    int idx_1;
    int idx_2;
    int myMax;
    int strt = 0;

    std::vector<int> complement;
    CompsDistinctSetup(z, complement, tar, idx_1, idx_2, nz, myMax);

    for (int i = width - nz, j = nz, nextStep = 0; i < width; ++i, --j) {

        nextStep += CountCompsDistinctLen(tar, i);

        CompsDistinctWorker(
            mat, v, z, complement, idx_1, idx_2,
            myMax, tar, strt, width, nextStep, nRows
        );

        strt = nextStep;

        std::iota(z.begin() + j - 1, z.end(), 1);
        z.back() = tar - static_cast<int>((i * (i + 1)) / 2);
        complement = PrepareComplement(z, tar);
    }

    CompsDistinctWorker(
        mat, v, z, complement, idx_1, idx_2,
        myMax, tar, strt, width, nRows, nRows
    );

    return 1;
}

template <typename T>
int CompsGenDistinct(
    RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
    std::vector<int> &z, std::size_t strt,
    std::size_t width, std::size_t nRows
) {

    int nz;
    int tar;
    int idx_1;
    int idx_2;
    int myMax;

    std::vector<int> complement;
    CompsDistinctSetup(z, complement, tar, idx_1, idx_2, nz, myMax);

    for (int i = width - nz, j = nz, nextStep = 0; i < width; ++i, --j) {

        nextStep += CountCompsDistinctLen(tar, i);

        CompsDistinctWorker(
            mat, v, z, complement, idx_1, idx_2,
            myMax, tar, strt, width, nextStep
        );

        strt = nextStep;

        std::iota(z.begin() + j - 1, z.end(), 1);
        z.back() = tar - static_cast<int>((i * (i + 1)) / 2);
        complement = PrepareComplement(z, tar);
    }

    CompsDistinctWorker(
        mat, v, z, complement, idx_1, idx_2,
        myMax, tar, strt, width, nRows
    );

    return 1;
}

int CompsDistinct(int* mat, std::vector<int> &z,
                  std::size_t width, std::size_t nRows) {

    int nz;
    int tar;
    int idx_1;
    int idx_2;
    int myMax;
    int strt = 0;

    std::vector<int> complement;
    CompsDistinctSetup(z, complement, tar, idx_1, idx_2, nz, myMax);

    for (int i = width - nz, j = nz, nextStep = 0; i < width; ++i, --j) {

        nextStep += CountCompsDistinctLen(tar, i);

        CompsDistinctWorker(
            mat, z, complement, idx_1, idx_2,
            myMax, tar, strt, width, nextStep, nRows
        );

        strt = nextStep;

        std::iota(z.begin() + j - 1, z.end(), 1);
        z.back() = tar - static_cast<int>((i * (i + 1)) / 2);
        complement = PrepareComplement(z, tar);
    }

    CompsDistinctWorker(
        mat, z, complement, idx_1, idx_2,
        myMax, tar, strt, width, nRows, nRows
    );

    return 1;
}

int CompsDistinct(
    RcppParallel::RMatrix<int> &mat, std::vector<int> &z,
    std::size_t strt, std::size_t width, std::size_t nRows
) {

    int nz;
    int tar;
    int idx_1;
    int idx_2;
    int myMax;

    std::vector<int> complement;
    CompsDistinctSetup(z, complement, tar, idx_1, idx_2, nz, myMax);

    for (int i = width - nz, j = nz, nextStep = 0; i < width; ++i, --j) {

        nextStep += CountCompsDistinctLen(tar, i);

        CompsDistinctWorker(
            mat, z, complement, idx_1, idx_2,
            myMax, tar, strt, width, nextStep
        );

        strt = nextStep;

        std::iota(z.begin() + j - 1, z.end(), 1);
        z.back() = tar - static_cast<int>((i * (i + 1)) / 2);
        complement = PrepareComplement(z, tar);
    }

    CompsDistinctWorker(
        mat, z, complement, idx_1, idx_2,
        myMax, tar, strt, width, nRows
    );

    return 1;
}

template int CompsGenDistinct(int*, const std::vector<int>&,
                              std::vector<int>&, std::size_t, std::size_t);
template int CompsGenDistinct(double*, const std::vector<double>&,
                              std::vector<int>&, std::size_t, std::size_t);

template int CompsGenDistinct(
    RcppParallel::RMatrix<int>&, const std::vector<int>&,
    std::vector<int>&, std::size_t, std::size_t, std::size_t
);

template int CompsGenDistinct(
    RcppParallel::RMatrix<double>&, const std::vector<double>&,
    std::vector<int>&, std::size_t, std::size_t, std::size_t
);
