#include "Partitions/CompositionsDistinctUtils.h"
#include "Partitions/PartitionsCountDistinct.h"
#include "Partitions/NextComposition.h"
#include <algorithm>
#include "RMatrix.h"
#include <numeric>
#include <limits>

constexpr int max_int = std::numeric_limits<int>::max();

// If the size of complement is 0 then we know we are simply taking
// permutations of z. Same thing with size 1 but a little reasoning is
// needed. In NextCompositionDistinct, we are swapping the last two
// elements in z with elements in complement. If complement only has one
// element, we are reduced to the case of taking permutations of z.
void PermutationsOnlyBranch(
    int* mat, std::vector<int> &z, std::size_t strt,
    std::size_t width, std::size_t rowLimit, std::size_t nRows
) {

    std::size_t count = strt;

    do {
        for (std::size_t k = 0; k < width; ++k) {
            mat[count + nRows * k] = z[k];
        }

        ++count;
    } while (std::next_permutation(z.begin(), z.end()) && count < rowLimit);
}

template <typename T>
void PermutationsOnlyBranch(
    T* mat, const std::vector<T> &v, std::vector<int> &z, std::size_t strt,
    std::size_t width, std::size_t rowLimit, std::size_t nRows
) {

    std::size_t count = strt;

    do {
        for (std::size_t k = 0; k < width; ++k) {
            mat[count + nRows * k] = v[z[k]];
        }

        ++count;
    } while (std::next_permutation(z.begin(), z.end()) && count < rowLimit);
}

void CompsDistStdWorker(
    int* mat, std::vector<int> &z, std::vector<int> &complement,
    int i1, int i2, int myMax, int tar, std::size_t strt,
    std::size_t width, std::size_t nRows
) {

    if (complement.size() < 2) {
        PermutationsOnlyBranch(mat, z, strt, width, nRows, nRows);
        return;
    }

    std::vector<int> idx;
    std::vector<int> tailSum;

    for (std::size_t count = strt, m = width - 1, q = complement.size() - 1,
         lastRow = nRows - 1; count < lastRow; ++count) {

        for (std::size_t k = 0; k < width; ++k) {
            mat[count + nRows * k] = z[k];
        }

        NextCompositionDistinct(
            z, complement, idx, tailSum, i1, i2, myMax, m, q, tar
        );
    }

    for (std::size_t k = 0; k < width; ++k) {
        mat[nRows - 1 + nRows * k] = z[k];
    }
}

template <typename T>
void CompsDistStdWorker(
    T* mat, const std::vector<T> &v, std::vector<int> &z,
    std::vector<int> &complement, int i1, int i2, int myMax,
    int tar, std::size_t strt, std::size_t width, std::size_t nRows
) {

    if (complement.size() < 2) {
        PermutationsOnlyBranch(mat, v, z, strt, width, nRows, nRows);
        return;
    }

    std::vector<int> idx;
    std::vector<int> tailSum;

    for (std::size_t count = strt, m = width - 1, q = complement.size() - 1,
         lastRow = nRows - 1; count < lastRow; ++count) {

        for (std::size_t k = 0; k < width; ++k) {
            mat[count + nRows * k] = v[z[k]];
        }

        NextCompositionDistinct(
            z, complement, idx, tailSum, i1, i2, myMax, m, q, tar
        );
    }

    for (std::size_t k = 0; k < width; ++k) {
        mat[nRows - 1 + nRows * k] = v[z[k]];
    }
}

void CompsDistStdWorker(
    RcppParallel::RMatrix<int> &mat, std::vector<int> &z,
    std::vector<int> &complement, int i1, int i2, int myMax, int tar,
    std::size_t strt, std::size_t width, std::size_t nRows
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
void CompsDistStdWorker(
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

void CompsDistMZWorker(
    int* mat, std::vector<int> &z, std::vector<int> &complement,
    int i1, int i2, int myMax, int tar, int nz, std::size_t strt,
    std::size_t width, int &rowLimit, std::size_t nRows
) {

    if (complement.size() < 2) {
        PermutationsOnlyBranch(mat, z, strt, width, rowLimit, nRows);
        return;
    }

    std::vector<int> idx;
    std::vector<int> tailSum;
    std::size_t count = strt;

    for (std::size_t m = width - 1, q = complement.size() - 1,
         lastRow = rowLimit - 1; count < lastRow && z[1]; ++count) {

        for (std::size_t k = 1; k < width; ++k) {
            mat[count + nRows * (k + nz - 1)] = z[k];
        }

        NextCompositionDistinct(
            z, complement, idx, tailSum, i1, i2, myMax, m, q, tar
        );
    }

    if (z[1]) {
        for (std::size_t k = 1; k < width && z[1]; ++k) {
            mat[rowLimit - 1 + nRows * (k + nz - 1)] = z[k];
        }
    } else {
        rowLimit = count;
    }
}

template <typename T>
void CompsDistMZWorker(
    T* mat, const std::vector<T> &v, std::vector<int> &z,
    std::vector<int> &complement, int i1, int i2, int myMax, int tar, int nz,
    std::size_t strt, std::size_t width, int &rowLimit, std::size_t nRows
) {

    if (complement.size() < 2) {
        PermutationsOnlyBranch(mat, v, z, strt, width, rowLimit, nRows);
        return;
    }

    std::vector<int> idx;
    std::vector<int> tailSum;
    std::size_t count = strt;

    for (std::size_t m = width - 1, q = complement.size() - 1,
         lastRow = rowLimit - 1; count < lastRow && z[1]; ++count) {

        for (std::size_t k = 1; k < width; ++k) {
            mat[count + nRows * (k + nz - 1)] = v[z[k]];
        }

        NextCompositionDistinct(
            z, complement, idx, tailSum, i1, i2, myMax, m, q, tar
        );
    }

    if (z[1]) {
        for (std::size_t k = 1; k < width && z[1]; ++k) {
            mat[rowLimit - 1 + nRows * (k + nz - 1)] = v[z[k]];
        }
    } else {
        rowLimit = count;
    }
}

void CompsDistMZWorker(
    RcppParallel::RMatrix<int> &mat, std::vector<int> &z,
    std::vector<int> &complement, int i1, int i2, int myMax, int tar,
    int nz, std::size_t strt, std::size_t width, int &nRows
) {

    std::vector<int> idx;
    std::vector<int> tailSum;
    std::size_t count = strt;

    for (std::size_t m = width - 1, q = complement.size() - 1,
         lastRow = nRows - 1; count < lastRow && z[1]; ++count) {

        for (std::size_t k = 1; k < width; ++k) {
            mat(count, k + nz - 1) = z[k];
        }

        NextCompositionDistinct(
            z, complement, idx, tailSum, i1, i2, myMax, m, q, tar
        );
    }

    if (z[1]) {
        for (std::size_t k = 1; k < width && z[1]; ++k) {
            mat(nRows - 1, k + nz - 1) = z[k];
        }
    } else {
        nRows = count;
    }
}

template <typename T>
void CompsDistMZWorker(
    RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
    std::vector<int> &z, std::vector<int> &complement,
    int i1, int i2, int myMax, int tar, int nz, std::size_t strt,
    std::size_t width, int &nRows
) {

    std::vector<int> idx;
    std::vector<int> tailSum;
    std::size_t count = strt;

    for (std::size_t m = width - 1, q = complement.size() - 1,
         lastRow = nRows - 1; count < lastRow && z[1]; ++count) {

        for (std::size_t k = 1; k < width; ++k) {
            mat(count, k + nz - 1) = v[z[k]];
        }

        NextCompositionDistinct(
            z, complement, idx, tailSum, i1, i2, myMax, m, q, tar
        );
    }

    if (z[1]) {
        for (std::size_t k = 1; k < width && z[1]; ++k) {
            mat(nRows - 1, k + nz - 1) = v[z[k]];
        }
    } else {
        nRows = count;
    }
}

int CompsDistinct(int* mat, std::vector<int> &z, std::size_t width,
                  std::size_t nRows, bool isWeak, int zeroBudget) {

    int tar = 0;
    int idx_1 = 0;
    int idx_2 = 0;
    int myMax = 0;

    std::size_t strt = 0;
    std::vector<int> complement;
    const int nz = std::count(z.cbegin(), z.cend(), 0);

    if (!isWeak && nz) {
        if (nz > 1) z.erase(z.begin(), z.begin() + (nz - 1));

        for (int i = width - nz, j = nz, nextStep = 0; i < width; ++i, --j) {

            CompsDistinctSetup(
                z, complement, tar, idx_1, idx_2, myMax, max_int, false, 1
            );

            double safe_count = CountCompsDistinctLen(tar, i);

            if ((safe_count + nextStep) < max_int) {
                nextStep += safe_count;
            } else {
                nextStep = max_int;
            }

            nextStep = std::min(nextStep, static_cast<int>(nRows));

            CompsDistMZWorker(
                mat, z, complement, idx_1, idx_2, myMax,
                tar, j, strt, z.size(), nextStep, nRows
            );

            for (std::size_t k = 0; k < j; ++k) {
                for (std::size_t count = strt,
                     offSet = nRows * k; count < nextStep; ++count) {
                    mat[count + offSet] = 0;
                }
            }

            if (nextStep >= static_cast<int>(nRows)) return 1;
            strt = nextStep;

            std::iota(z.begin(), z.end(), 1);
            z.back() = tar - static_cast<int>((i * (i + 1)) / 2);
            if (j > 1) z.insert(z.begin(), 0);
        }
    }

    zeroBudget = isWeak ? zeroBudget : std::count(z.cbegin(), z.cend(), 0);
    CompsDistinctSetup(
        z, complement, tar, idx_1, idx_2, myMax, max_int, isWeak, zeroBudget
    );
    CompsDistStdWorker(
        mat, z, complement, idx_1, idx_2, myMax, tar, strt, width, nRows
    );

    return 1;
}

template <typename T>
int CompsGenDistinct(
    T* mat, const std::vector<T> &v, std::vector<int> &z,
    std::size_t width, std::size_t nRows, bool isWeak, int zeroBudget
) {

    int tar = 0;
    int idx_1 = 0;
    int idx_2 = 0;
    int myMax = 0;

    std::size_t strt = 0;
    std::vector<int> complement;
    const int nz = std::count(z.cbegin(), z.cend(), 0);
    const int idx_max = static_cast<int>(v.size()) - 1;

    if (!isWeak && nz && v.front() == 0) {
        if (nz > 1) z.erase(z.begin(), z.begin() + (nz - 1));
        std::vector<int> allowed(idx_max);
        std::iota(allowed.begin(), allowed.end(), 1);

        for (int i = width - nz, j = nz, nextStep = 0; i < width; ++i, --j) {

            CompsDistinctSetup(
                z, complement, tar, idx_1, idx_2, myMax, idx_max, false, 1
            );

            double safe_count = CountCompDistLenRstrctd(tar, i, allowed);

            if ((safe_count + nextStep) < max_int) {
                nextStep += safe_count;
            } else {
                nextStep = max_int;
            }

            nextStep = std::min(nextStep, static_cast<int>(nRows));

            CompsDistMZWorker(
                mat, v, z, complement, idx_1, idx_2,
                myMax, tar, j, strt, z.size(), nextStep, nRows
            );

            for (std::size_t k = 0; k < j; ++k) {
                for (std::size_t count = strt,
                     offSet = nRows * k; count < nextStep; ++count) {
                    mat[count + offSet] = v[0];
                }
            }

            if (nextStep >= static_cast<int>(nRows)) return 1;
            strt = nextStep;

            if (idx_max == tar) {
                std::iota(z.begin(), z.end(), 1);
                z.back() = tar - static_cast<int>((i * (i + 1)) / 2);
            } else {
                GetFirstPartitionDistinct(allowed, z, tar, z.size(), idx_max);
                for (auto &z_i: z) ++z_i;
            }

            if (j > 1) z.insert(z.begin(), 0);
        }
    }

    bool compZero = IsComplementZeroBased(v.front() == 0, isWeak, true);
    zeroBudget = isWeak ? zeroBudget : std::count(z.cbegin(), z.cend(), 0);

    CompsDistinctSetup(
        z, complement, tar, idx_1, idx_2, myMax, idx_max, compZero, zeroBudget
    );
    CompsDistStdWorker(
        mat, v, z, complement, idx_1, idx_2, myMax, tar, strt, width, nRows
    );

    return 1;
}

int CompsDistinct(
    RcppParallel::RMatrix<int> &mat, std::vector<int> &z,
    std::size_t strt, std::size_t width, std::size_t nRows,
    bool isWeak, int zeroBudget
) {

    int tar = 0;
    int idx_1 = 0;
    int idx_2 = 0;
    int myMax = 0;

    std::vector<int> complement;
    const int nz = std::count(z.cbegin(), z.cend(), 0);

    if (!isWeak && nz) {
        if (nz > 1) z.erase(z.begin(), z.begin() + (nz - 1));

        for (int i = width - nz, j = nz, nextStep = 0; i < width; ++i, --j) {

            CompsDistinctSetup(
                z, complement, tar, idx_1, idx_2, myMax, max_int, false, 1
            );

            double safe_count = CountCompsDistinctLen(tar, i);

            if ((safe_count + nextStep) < max_int) {
                nextStep += safe_count;
            } else {
                nextStep = max_int;
            }

            nextStep = std::min(nextStep, static_cast<int>(nRows));

            CompsDistMZWorker(
                mat, z, complement, idx_1, idx_2,
                myMax, tar, j, strt, z.size(), nextStep
            );

            for (std::size_t k = 0; k < j; ++k) {
                for (std::size_t count = strt; count < nextStep; ++count) {
                    mat(count, k) = 0;
                }
            }

            if (nextStep >= static_cast<int>(nRows)) return 1;
            strt = nextStep;

            std::iota(z.begin(), z.end(), 1);
            z.back() = tar - static_cast<int>((i * (i + 1)) / 2);
            if (j > 1) z.insert(z.begin(), 0);
        }
    }

    zeroBudget = isWeak ? zeroBudget : std::count(z.cbegin(), z.cend(), 0);
    CompsDistinctSetup(
        z, complement, tar, idx_1, idx_2, myMax, max_int, isWeak, zeroBudget
    );
    CompsDistStdWorker(
        mat, z, complement, idx_1, idx_2, myMax, tar, strt, width, nRows
    );

    return 1;
}

template <typename T>
int CompsGenDistinct(
    RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
    std::vector<int> &z, std::size_t strt, std::size_t width,
    std::size_t nRows, bool isWeak, int zeroBudget
) {

    int tar = 0;
    int idx_1 = 0;
    int idx_2 = 0;
    int myMax = 0;

    std::vector<int> complement;
    const int nz = std::count(z.cbegin(), z.cend(), 0);
    const int idx_max = static_cast<int>(v.size()) - 1;

    if (!isWeak && nz && v.front() == 0) {
        if (nz > 1) z.erase(z.begin(), z.begin() + (nz - 1));
        std::vector<int> allowed(idx_max);
        std::iota(allowed.begin(), allowed.end(), 1);

        for (int i = width - nz, j = nz, nextStep = 0; i < width; ++i, --j) {

            CompsDistinctSetup(
                z, complement, tar, idx_1, idx_2, myMax, idx_max, false, 1
            );

            double safe_count = CountCompDistLenRstrctd(tar, i, allowed);

            if ((safe_count + nextStep) < max_int) {
                nextStep += safe_count;
            } else {
                nextStep = max_int;
            }

            nextStep = std::min(nextStep, static_cast<int>(nRows));

            CompsDistMZWorker(
                mat, v, z, complement, idx_1, idx_2,
                myMax, tar, j, strt, z.size(), nextStep
            );

            for (std::size_t k = 0; k < j; ++k) {
                for (std::size_t count = strt; count < nextStep; ++count) {
                    mat(count, k) = v[0];
                }
            }

            if (nextStep >= static_cast<int>(nRows)) return 1;
            strt = nextStep;

            if (idx_max == tar) {
                std::iota(z.begin(), z.end(), 1);
                z.back() = tar - static_cast<int>((i * (i + 1)) / 2);
            } else {
                GetFirstPartitionDistinct(allowed, z, tar, z.size(), idx_max);
                for (auto &z_i: z) ++z_i;
            }

            if (j > 1) z.insert(z.begin(), 0);
        }
    }

    bool compZero = IsComplementZeroBased(v.front() == 0, isWeak, true);
    zeroBudget = isWeak ? zeroBudget : std::count(z.cbegin(), z.cend(), 0);

    CompsDistinctSetup(
        z, complement, tar, idx_1, idx_2, myMax, idx_max, compZero, zeroBudget
    );
    CompsDistStdWorker(
        mat, v, z, complement, idx_1, idx_2, myMax, tar, strt, width, nRows
    );

    return 1;
}

template int CompsGenDistinct(int*, const std::vector<int>&, std::vector<int>&,
                              std::size_t, std::size_t, bool, int);
template int CompsGenDistinct(
    double*, const std::vector<double>&, std::vector<int>&,
    std::size_t, std::size_t, bool, int
);

template int CompsGenDistinct(
    RcppParallel::RMatrix<int>&, const std::vector<int>&,
    std::vector<int>&, std::size_t, std::size_t, std::size_t, bool, int
);

template int CompsGenDistinct(
    RcppParallel::RMatrix<double>&, const std::vector<double>&,
    std::vector<int>&, std::size_t, std::size_t, std::size_t, bool, int
);
