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
void PartsGenMultiset(std::vector<T> &partitionsVec, const std::vector<T> &v,
                      const std::vector<int> &Reps, std::vector<int> &z,
                      int width, int lastElem, int lastCol, int maxRows) {

    int b = 0;
    int p = 0;
    int e = 0;

    std::vector<int> rpsCnt;
    PrepareMultisetPart(z, rpsCnt, b, p, e, lastCol, lastElem);

    for (int count = 0; keepGoing(rpsCnt, lastElem, z, e, b);
         NextMultisetGenPart(rpsCnt, z, e, b,  p, lastCol, lastElem)) {
        
        for (int k = 0; k < width; ++k) {
            partitionsVec.push_back(v[z[k]]);
        }

        ++count;

        if (count >= maxRows) {break;}
    }

    const int numResult = partitionsVec.size() / width;

    if (numResult < maxRows) {
        for (int k = 0; k < width; ++k) {
            partitionsVec.push_back(v[z[k]]);
        }
    }
}

template <typename T>
void PartsGenPermMultiset(std::vector<T> &partitionsVec,
                          const std::vector<T> &v,
                          const std::vector<int> &Reps, std::vector<int> &z,
                          int width, int lastElem, int lastCol, int maxRows) {

    int b = 0;
    int p = 0;
    int e = 0;

    std::vector<int> rpsCnt;
    PrepareMultisetPart(z, rpsCnt, b, p, e, lastCol, lastElem);

    for (int count = 0; keepGoing(rpsCnt, lastElem, z, e, b);
         NextMultisetGenPart(rpsCnt, z, e, b,  p, lastCol, lastElem)) {
        
        PopulateVecPerm(v, partitionsVec, z, count, width, maxRows);
        if (count >= maxRows) {break;}
    }

    int count = partitionsVec.size() / width;

    if (count < maxRows) {
        PopulateVecPerm(v, partitionsVec, z, count, width, maxRows);
    }
}

template void PartsGenMultiset(std::vector<int>&, const std::vector<int>&,
                               const std::vector<int>&, std::vector<int>&,
                               int, int, int, int);

template void PartsGenMultiset(std::vector<double>&,
                               const std::vector<double>&,
                               const std::vector<int>&, std::vector<int>&,
                               int, int, int, int);

template void PartsGenPermMultiset(std::vector<int>&, const std::vector<int>&,
                                   const std::vector<int>&, std::vector<int>&,
                                   int, int, int, int);

template void PartsGenPermMultiset(std::vector<double>&,
                                   const std::vector<double>&,
                                   const std::vector<int>&, std::vector<int>&,
                                   int, int, int, int);
