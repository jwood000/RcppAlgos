#include "Partitions/PopulateVecPerm.h"
#include "Partitions/NextPartition.h"

bool keepGoing(const std::vector<int> &rpsCnt, int lastElem,
               const std::vector<int> &z, int edge, int boundary) {

    if (edge >= 0) {
        const int myDiff = z[boundary] - z[edge];

        if (myDiff < 2) {
            return false;
        } else if (myDiff == 2) {
            return (rpsCnt[z[edge] + 1] > 1);
        } else {
            return (rpsCnt[z[edge] + 1] && rpsCnt[z[boundary] - 1]);
        }
    } else {
        return false;
    }
}

template <typename T>
void PartsGenMultiset(std::vector<T> &partsVec, const std::vector<T> &v,
                      const std::vector<int> &Reps, std::vector<int> &z,
                      int width, int nRows) {

    int b = 0;
    int p = 0;
    int e = 0;

    const int lastCol = width - 1;
    const int lastElem = v.size() - 1;
    std::vector<int> rpsCnt(Reps.cbegin(), Reps.cend());

    PrepareMultisetPart(rpsCnt, z, b, p, e, lastCol, lastElem);

    for (int count = 0; keepGoing(rpsCnt, lastElem, z, e, b);
         NextMultisetGenPart(rpsCnt, z, e, b, p, lastCol, lastElem)) {

        for (int k = 0; k < width; ++k) {
            partsVec.push_back(v[z[k]]);
        }

        ++count;

        if (count >= nRows) {break;}
    }

    const int numResult = partsVec.size() / width;

    if (numResult < nRows) {
        for (int k = 0; k < width; ++k) {
            partsVec.push_back(v[z[k]]);
        }
    }
}

template <typename T>
void PartsGenPermMultiset(std::vector<T> &partsVec,
                          const std::vector<T> &v,
                          const std::vector<int> &Reps,
                          std::vector<int> &z, int width, int nRows) {

    int b = 0;
    int p = 0;
    int e = 0;

    const int lastCol = width - 1;
    const int lastElem = v.size() - 1;
    std::vector<int> rpsCnt(Reps.cbegin(), Reps.cend());

    PrepareMultisetPart(rpsCnt, z, b, p, e, lastCol, lastElem);

    for (int count = 0; keepGoing(rpsCnt, lastElem, z, e, b);
         NextMultisetGenPart(rpsCnt, z, e, b,  p, lastCol, lastElem)) {

        PopulateVecPerm(v, partsVec, z, count, width, nRows);
        if (count >= nRows) {break;}
    }

    int count = partsVec.size() / width;

    if (count < nRows) {
        PopulateVecPerm(v, partsVec, z, count, width, nRows);
    }
}

template void PartsGenMultiset(std::vector<int>&, const std::vector<int>&,
                               const std::vector<int>&, std::vector<int>&,
                               int, int);

template void PartsGenMultiset(std::vector<double>&,
                               const std::vector<double>&,
                               const std::vector<int>&, std::vector<int>&,
                               int, int);

template void PartsGenPermMultiset(std::vector<int>&, const std::vector<int>&,
                                   const std::vector<int>&, std::vector<int>&,
                                   int, int);

template void PartsGenPermMultiset(std::vector<double>&,
                                   const std::vector<double>&,
                                   const std::vector<int>&, std::vector<int>&,
                                   int, int);
