#include "Partitions/NextPartition.h"
#include "PopulateVec.h"

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
                      std::size_t width, std::size_t nRows, bool IsComb) {

    int b = 0;
    int p = 0;
    int e = 0;

    const int lastCol  = width - 1;
    const int lastElem = v.size() - 1;
    std::vector<int> rpsCnt(Reps.cbegin(), Reps.cend());
    PrepareMultisetPart(rpsCnt, z, b, p, e, lastCol, lastElem);

    for (std::size_t count = 0; keepGoing(rpsCnt, lastElem, z, e, b);
         NextMultisetGenPart(rpsCnt, z, e, b, p, lastCol, lastElem)) {

        PopulateVec(v, partsVec, z, count, width, nRows, IsComb);
        if (count >= nRows) break;
    }

    std::size_t count = partsVec.size() / width;

    if (count < nRows) {
        PopulateVec(v, partsVec, z, count, width, nRows, IsComb);
    }
}

template void PartsGenMultiset(std::vector<int>&, const std::vector<int>&,
                               const std::vector<int>&, std::vector<int>&,
                               std::size_t, std::size_t, bool);
template void PartsGenMultiset(std::vector<double>&,
                               const std::vector<double>&,
                               const std::vector<int>&, std::vector<int>&,
                               std::size_t, std::size_t, bool);
