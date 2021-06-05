#include "Partitions/PopulateVecPerm.h"
#include "Partitions/NextPartition.h"

template <typename T>
void PartsGenRep(T* mat, const std::vector<T> &v, std::vector<int> &z,
                 int m, int lastElem, int lastCol, int strt, int nRows) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;

    PrepareRepPart(z, boundary, pivot, edge, lastElem, lastCol);

    for (int count = strt; count < nRows; ++count) {
        for (int k = 0; k < m; ++k) {
            mat[count + nRows * k] = v[z[k]];
        }
        
        NextRepGenPart(z, boundary, edge, pivot, lastCol, lastElem);
    }
}

template <typename T>
void PartsGenPermRep(std::vector<T> &partitionsVec,
                     const std::vector<T> &v, std::vector<int> &z,
                     int m, int lastElem, int lastCol, int maxRows) {

    int edge = 0;
    int pivot = 0;
    int count = 0;
    int boundary = 0;
    PrepareRepPart(z, boundary, pivot, edge, lastElem, lastCol);

    while ((edge >= 0) && (z[boundary] - z[edge] >= 2)) {
        PopulateVecPerm(v, partitionsVec, z, count, m, maxRows);

        if (count >= maxRows)
            break;

        NextRepGenPart(z, boundary, edge, pivot, lastCol, lastElem);
    }

    if (count < maxRows)
        PopulateVecPerm(v, partitionsVec, z, count, m, maxRows);
}

void PartsRep(int* mat, std::vector<int> &z, int m, int boundary,
              int edge, int lastCol, int strt, int nRows) {

    for (int count = strt; count < nRows; ++count,
         NextRepPart(z, boundary, edge, lastCol)) {

        for (std::size_t k = 0; k < m; ++k) {
            mat[count + nRows * k] = z[k];
        }
    }
}

void PartsPermRep(int* mat, std::vector<int> &z, int m,
                  int boundary, int edge, int lastCol, int nRows) {

    for (int count = 0; ;
         NextRepPart(z, boundary, edge, lastCol)) {

        do {
            for (int k = 0; k < m; ++k) {
                mat[count + nRows * k] = z[k];
            }
            
            ++count;
        } while (std::next_permutation(z.begin(), z.end()) && count < nRows);

        if (count >= nRows) {break;}
    }
}

template void PartsGenRep(int*, const std::vector<int>&,
                          std::vector<int>&, int, int, int, int, int);

template void PartsGenRep(double*, const std::vector<double>&,
                          std::vector<int>&, int, int, int, int, int);

template void PartsGenPermRep(std::vector<int>&, const std::vector<int>&,
                              std::vector<int>&, int, int, int, int);

template void PartsGenPermRep(std::vector<double>&, const std::vector<double>&,
                              std::vector<int>&, int, int, int, int);
